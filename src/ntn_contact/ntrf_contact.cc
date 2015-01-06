/**
 * @file   ntrf_contact.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// simtools
#include "ntrf_contact.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTRFContact::NTRFContact(SolidMechanicsModel & model,
			 const ContactID & id,
			 const MemoryID & memory_id) :
  NTNBaseContact(model,id,memory_id),
  reference_point(model.getSpatialDimension()),
  normal(model.getSpatialDimension())
{
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::setReferencePoint(Real x, Real y, Real z) {
  AKANTU_DEBUG_IN();

  Real coord[3];
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;

  UInt dim = this->model.getSpatialDimension();
  for (UInt d=0; d<dim; ++d)
    this->reference_point(d) = coord[d];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::setNormal(Real x, Real y, Real z) {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();

  Real coord[3];
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;

  for (UInt d=0; d<dim; ++d)
    this->normal(d) = coord[d];
  
  this->normal.normalize();
  
  this->updateNormals();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addSurface(const Surface & surf) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh_ref = this->model.getFEEngine().getMesh();

  try {
    const ElementGroup & boundary = mesh_ref.getElementGroup(surf);
    this->contact_surfaces.insert(&boundary);

    // find slave nodes
    for(ElementGroup::const_node_iterator nodes_it(boundary.node_begin());
	nodes_it!= boundary.node_end();
	++nodes_it) {
      if (!this->model.isPBCSlaveNode(*nodes_it)) {
	this->addSplitNode(*nodes_it);
      }
    }
  } catch (debug::Exception e) {
    AKANTU_DEBUG_INFO("NTRFContact addSurface did not found subboundary " << surf 
		      << " and ignored it. Other procs might have it :)");
  }

  // synchronize with depending nodes
  findBoundaryElements(this->slaves.getArray(), this->slave_elements);
  updateInternalData();
  syncArrays(_added);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addNodes(Array<UInt> & nodes) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = nodes.getSize();
  UInt nb_compo = nodes.getNbComponent();
  for (UInt n=0; n<nb_nodes; ++n) {
    for (UInt c=0; c<nb_compo; ++c) {
      this->addSplitNode(nodes(n,c));
    }
  }

  // synchronize with depending nodes
  findBoundaryElements(this->slaves.getArray(), this->slave_elements);
  updateInternalData();
  syncArrays(_added);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::updateNormals() {
  AKANTU_DEBUG_IN();
  
  // normal is the same for all slaves
  this->normals.set(this->normal);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::updateImpedance() {
  AKANTU_DEBUG_IN();
  
  UInt nb_contact_nodes = getNbContactNodes();
  Real delta_t = this->model.getTimeStep();
  AKANTU_DEBUG_ASSERT(delta_t != NAN, "Time step is NAN. Have you set it already?");

  const Array<Real> & mass = this->model.getMass();

  for (UInt n=0; n<nb_contact_nodes; ++n) {  
    UInt slave = this->slaves(n);

    Real imp = this->lumped_boundary_slaves(n) / mass(slave);
    imp = 2 / delta_t / imp;
    this->impedance(n) = imp;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeRelativeTangentialField(const Array<Real> & field,
						 Array<Real> & rel_tang_field) const {
  AKANTU_DEBUG_IN();

  // resize arrays to zero
  rel_tang_field.resize(0);

  UInt dim = this->model.getSpatialDimension();
  
  Array<Real>::const_iterator< Vector<Real> > it_field  = field.begin(dim);
  Array<Real>::const_iterator< Vector<Real> > it_normal = this->normals.getArray().begin(dim);

  Vector<Real> rfv(dim);
  Vector<Real> np_rfv(dim);

  UInt nb_contact_nodes = this->slaves.getSize();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt node = this->slaves(n);

    // relative field vector
    rfv = it_field[node];;

    // normal projection of relative field
    const Vector<Real> & normal_v = it_normal[n];
    np_rfv = normal_v;
    np_rfv *= rfv.dot(normal_v);

    // subtract normal projection from relative field to get the tangential projection
    rfv -= np_rfv;
    rel_tang_field.push_back(rfv);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeNormalGap(Array<Real> & gap) const {
  AKANTU_DEBUG_IN();

  gap.resize(0);

  UInt dim = this->model.getSpatialDimension();

  Array<Real>::const_iterator< Vector<Real> > it_cur_pos = this->model.getCurrentPosition().begin(dim);
  Array<Real>::const_iterator< Vector<Real> > it_normal = this->normals.getArray().begin(dim);

  Vector<Real> gap_v(dim);

  UInt nb_contact_nodes = this->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt node = this->slaves(n);

    // gap vector
    gap_v = it_cur_pos[node];
    gap_v -= this->reference_point;

    // length of normal projection of gap vector
    const Vector<Real> & normal_v  = it_normal[n];
    gap.push_back(gap_v.dot(normal_v));
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeRelativeNormalField(const Array<Real> & field,
					     Array<Real> & rel_normal_field) const {
  AKANTU_DEBUG_IN();

  // resize arrays to zero
  rel_normal_field.resize(0);

  UInt dim = this->model.getSpatialDimension();

  Array<Real>::const_iterator< Vector<Real> > it_field  = field.begin(dim);
  Array<Real>::const_iterator< Vector<Real> > it_normal = this->normals.getArray().begin(dim);

  UInt nb_contact_nodes = this->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt node = this->slaves(n);

    const Vector<Real> & field_v  = it_field[node];
    const Vector<Real> & normal_v = it_normal[n];
    
    rel_normal_field.push_back(field_v.dot(normal_v));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "NTRFContact [" << std::endl;
  NTNBaseContact::printself(stream,indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addDumpFieldToDumper(const std::string & dumper_name,
				       const std::string & field_id) {
  AKANTU_DEBUG_IN();

  /*
#ifdef AKANTU_USE_IOHELPER
  const Array<UInt> & nodal_filter = this->slaves.getArray();

#define ADD_FIELD(field_id, field, type)				\
  internalAddDumpFieldToDumper(dumper_name,				\
			       field_id,				\
			       new DumperIOHelper::NodalField< type, true, \
							       Array<type>, \
							       Array<UInt> >(field, 0, 0, &nodal_filter))
  */

  /*
  if(field_id == "displacement") {
    ADD_FIELD(field_id, this->model.getDisplacement(), Real);
  }
  else if(field_id == "contact_pressure") {
    internalAddDumpFieldToDumper(dumper_name,
				 field_id,
				 new DumperIOHelper::NodalField<Real>(this->contact_pressure.getArray()));
  }
  else {*/
  NTNBaseContact::addDumpFieldToDumper(dumper_name, field_id);
    //}

  /*
#undef ADD_FIELD
#endif
  */

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
