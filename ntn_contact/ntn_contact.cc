/**
 * @file   ntn_contact.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Mar 14 11:52:00 2013
 *
 * @brief  implementation of ntn_contact
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_contact.hh"
#include "dumper_text.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTNContact::NTNContact(SolidMechanicsModel & model,
                       const ContactID & id,
                       const MemoryID & memory_id) : 
  NTNBaseContact(model,id,memory_id),
  masters(0,1,0,id+":masters",std::numeric_limits<UInt>::quiet_NaN(),"masters"),
  lumped_boundary_masters(0,1,0,id+":lumped_boundary_masters",
			  std::numeric_limits<Real>::quiet_NaN(),"lumped_boundary_masters"),
  master_elements("master_elements", id, memory_id)
{
  AKANTU_DEBUG_IN();
  
  const Mesh & mesh = this->model.getMesh();
  UInt spatial_dimension = this->model.getSpatialDimension();
  
  mesh.initElementTypeMapArray(this->master_elements,
			       1,
			       spatial_dimension - 1);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::pairInterfaceNodes(const ElementGroup & slave_boundary, 
				    const ElementGroup & master_boundary,
				    UInt surface_normal_dir,
				    const Mesh & mesh,
				    Array<UInt> & pairs) {
  AKANTU_DEBUG_IN();

  pairs.resize(0);
  AKANTU_DEBUG_ASSERT(pairs.getNbComponent() == 2,
		      "Array of node pairs should have nb_component = 2," <<
		      " but has nb_component = " << pairs.getNbComponent());

  UInt dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(surface_normal_dir < dim, "Mesh is of " << dim << " dimensions"
		      << " and cannot have direction " << surface_normal_dir
		      << " for surface normal");
  
  // offset for projection computation
  UInt offset[dim-1];
  for (UInt i=0, j=0; i<dim; ++i) {
    if (surface_normal_dir != i) {
      offset[j] = i;
      ++j;
    }
  }

  // find projected node coordinates
  const Array<Real> & coordinates = mesh.getNodes();
  
  // find slave nodes
  Array<Real> proj_slave_coord(slave_boundary.getNbNodes(),dim-1,0.);
  Array<UInt> slave_nodes(slave_boundary.getNbNodes());
  UInt n(0);
  for(ElementGroup::const_node_iterator nodes_it(slave_boundary.node_begin()); nodes_it!= slave_boundary.node_end(); ++nodes_it) {
    for (UInt d=0; d<dim-1; ++d) {
      proj_slave_coord(n,d) = coordinates(*nodes_it,offset[d]);
      slave_nodes(n)=*nodes_it;
    }
    ++n;
  }
  
  // find master nodes
  Array<Real> proj_master_coord(master_boundary.getNbNodes(),dim-1,0.);
  Array<UInt> master_nodes(master_boundary.getNbNodes());
  n=0;
  for(ElementGroup::const_node_iterator nodes_it(master_boundary.node_begin()); nodes_it!= master_boundary.node_end(); ++nodes_it) {
    for (UInt d=0; d<dim-1; ++d) {
      proj_master_coord(n,d) = coordinates(*nodes_it,offset[d]);
      master_nodes(n)=*nodes_it;
    }
    ++n;
  }

  // find minimum distance between slave nodes to define tolerance
  Real min_dist = std::numeric_limits<Real>::max();
  for (UInt i=0; i<proj_slave_coord.getSize(); ++i) {
    for (UInt j=i+1; j<proj_slave_coord.getSize(); ++j) {
      Real dist = 0.;
      for (UInt d=0; d<dim-1; ++d) {
	dist += (proj_slave_coord(i,d) - proj_slave_coord(j,d)) 
	      * (proj_slave_coord(i,d) - proj_slave_coord(j,d));
      }
      if (dist < min_dist) {
	min_dist = dist;
      }
    }
  }
  min_dist = std::sqrt(min_dist);
  Real local_tol = 0.1*min_dist;

  // find master slave node pairs
  for (UInt i=0; i<proj_slave_coord.getSize(); ++i) {
    for (UInt j=0; j<proj_master_coord.getSize(); ++j) {
      Real dist = 0.;
      for (UInt d=0; d<dim-1; ++d) {
	dist += (proj_slave_coord(i,d) - proj_master_coord(j,d)) 
	      * (proj_slave_coord(i,d) - proj_master_coord(j,d));
      }
      dist = std::sqrt(dist);
      if (dist < local_tol) { // it is a pair
	UInt pair[2];
	pair[0] = slave_nodes(i);
	pair[1] = master_nodes(j);
	pairs.push_back(pair);
	continue; // found master do not need to search further for this slave
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addSurfacePair(const Surface & slave, 
				const Surface & master, 
				UInt surface_normal_dir) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = this->model.getMesh();

  const ElementGroup & slave_boundary  = mesh.getElementGroup(slave);
  const ElementGroup & master_boundary = mesh.getElementGroup(master);
  
  this->contact_surfaces.insert(&slave_boundary);
  this->contact_surfaces.insert(&master_boundary);

  Array<UInt> pairs(0,2);
  NTNContact::pairInterfaceNodes(slave_boundary, 
				 master_boundary, 
				 surface_normal_dir,
				 this->model.getMesh(),
				 pairs);
  
  this->addNodePairs(pairs);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addNodePairs(const Array<UInt> & pairs) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_ASSERT(pairs.getNbComponent() == 2, 
		      "Array of node pairs should have nb_component = 2," << 
		      " but has nb_component = " << pairs.getNbComponent());
  UInt nb_pairs = pairs.getSize();
  for (UInt n=0; n<nb_pairs; ++n) {
    this->addSplitNode(pairs(n,0), pairs(n,1));
  }
  
  // synchronize with depending nodes
  findBoundaryElements(this->slaves.getArray(),  this->slave_elements);
  findBoundaryElements(this->masters.getArray(), this->master_elements);
  updateInternalData();
  syncArrays(_added);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::getNodePairs(Array<UInt> & pairs) const {
  AKANTU_DEBUG_IN();
  
  pairs.resize(0);
  AKANTU_DEBUG_ASSERT(pairs.getNbComponent() == 2,
		      "Array of node pairs should have nb_component = 2," <<
		      " but has nb_component = " << pairs.getNbComponent());
  UInt nb_pairs = this->getNbContactNodes();
  for (UInt n=0; n<nb_pairs; ++n) {
    UInt pair[2];
    pair[0] = this->slaves(n);
    pair[1] = this->masters(n);
    pairs.push_back(pair);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addSplitNode(UInt slave, UInt master) {
  AKANTU_DEBUG_IN();

  NTNBaseContact::addSplitNode(slave);
  
  this->masters.push_back(master);
  this->lumped_boundary_masters.push_back(std::numeric_limits<Real>::quiet_NaN());
  
  AKANTU_DEBUG_OUT();  
}

/* -------------------------------------------------------------------------- */
/*
  This function only works for surface elements with one quad point. For 
  surface elements with more quad points, it computes still, but the result 
  might not be what you are looking for.
 */
void NTNContact::updateNormals() {
  AKANTU_DEBUG_IN();

  // set normals to zero
  this->normals.clear();

  // contact information
  UInt dim = this->model.getSpatialDimension();
  UInt nb_contact_nodes = this->getNbContactNodes();

  this->synch_registry->synchronize(_gst_cf_nodal); // synchronize current pos
  const Array<Real> & cur_pos = this->model.getCurrentPosition();

  FEEngine & boundary_fem = this->model.getFEEngineBoundary();
  const Mesh & mesh = this->model.getMesh();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    Mesh::type_iterator it   = mesh.firstType(dim-1, *gt);
    Mesh::type_iterator last = mesh.lastType(dim-1, *gt);

    for (; it != last; ++it) {
      // get elements connected to each node
      CSR<UInt> node_to_element;
      MeshUtils::buildNode2ElementsElementTypeMap(mesh, node_to_element, *it, *gt);
      
      // compute the normals
      Array<Real> quad_normals(0,dim);
      boundary_fem.computeNormalsOnControlPoints(cur_pos, quad_normals, *it, *gt);
      
      UInt nb_quad_points = boundary_fem.getNbQuadraturePoints(*it, *gt);
      
      // add up normals to all master nodes
      for (UInt n=0; n<nb_contact_nodes; ++n) {
	UInt master = this->masters(n);
	CSR<UInt>::iterator elem = node_to_element.begin(master);
	// loop over all elements connected to this master node
	for (; elem != node_to_element.end(master); ++elem) {
	  UInt e = *elem;
	  // loop over all quad points of this element
	  for (UInt q=0; q<nb_quad_points; ++q) {
	    // add quad normal to master normal
	    for (UInt d=0; d<dim; ++d) {
	      this->normals(n,d) += quad_normals(e*nb_quad_points + q, d);
	    }
	  }
	}
      }
    }
  }

  Real * master_normals = this->normals.storage();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    if (dim==2) 
      Math::normalize2(&(master_normals[n*dim]));
    else if (dim==3)
      Math::normalize3(&(master_normals[n*dim]));
  }

  // // normalize normals
  // Array<Real>::iterator<Real> nit  = this->normals.begin();
  // Array<Real>::iterator<Real> nend = this->normals.end();
  // for (; nit != nend; ++nit) {
  //   nit->normalize();
  // }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  NTNBaseContact::dumpRestart(file_name);
  this->masters.dumpRestartFile(file_name);
  this->lumped_boundary_masters.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  NTNBaseContact::readRestart(file_name);
  this->masters.readRestartFile(file_name);
  this->lumped_boundary_masters.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::updateImpedance() {
  AKANTU_DEBUG_IN();

  UInt nb_contact_nodes = getNbContactNodes();
  Real delta_t = this->model.getTimeStep();
  AKANTU_DEBUG_ASSERT(delta_t != NAN, "Time step is NAN. Have you set it already?");

  const Array<Real> & mass = this->model.getMass();

  for (UInt n=0; n<nb_contact_nodes; ++n) {  
    UInt master = this->masters(n);
    UInt slave = this->slaves(n);

    Real imp = (this->lumped_boundary_masters(n) / mass(master)) 
             + (this->lumped_boundary_slaves(n) / mass(slave));
    imp = 2 / delta_t / imp;
    this->impedance(n) = imp;
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::updateLumpedBoundary() {
  AKANTU_DEBUG_IN();

  internalUpdateLumpedBoundary(this->slaves.getArray(),
			       this->slave_elements,
			       this->lumped_boundary_slaves);
  
  internalUpdateLumpedBoundary(this->masters.getArray(),
			       this->master_elements,
			       this->lumped_boundary_masters);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::applyContactPressure() {
  AKANTU_DEBUG_IN();

  UInt nb_ntn_pairs = getNbContactNodes();
  UInt dim = this->model.getSpatialDimension();

  Array<Real> & residual = this->model.getResidual();

  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    UInt master = this->masters(n);
    UInt slave = this->slaves(n);
    
    for (UInt d=0; d<dim; ++d) {
      residual(master,d) += this->lumped_boundary_masters(n) * this->contact_pressure(n,d);
      residual(slave, d) -= this->lumped_boundary_slaves(n)  * this->contact_pressure(n,d);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::computeRelativeTangentialField(const Array<Real> & field,
						Array<Real> & rel_tang_field) const {
  AKANTU_DEBUG_IN();
  
  // resize arrays to zero
  rel_tang_field.resize(0);

  UInt dim = this->model.getSpatialDimension();

  Array<Real>::const_vector_iterator it_field  = field.begin(dim);
  Array<Real>::const_vector_iterator it_normal = this->normals.getArray().begin(dim);

  Vector<Real> rfv(dim);
  Vector<Real> np_rfv(dim);

  UInt nb_contact_nodes = this->slaves.getSize();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt slave  = this->slaves(n);
    UInt master = this->masters(n);

    // relative field vector (slave - master)
    rfv = it_field[slave];
    rfv -= it_field[master];

    // normal projection of relative field
    const Vector<Real> & normal_v = it_normal[n];
    np_rfv = normal_v;
    np_rfv *= rfv.dot(normal_v);

    // subract normal projection from relative field to get the tangential projection
    rfv -= np_rfv;
    rel_tang_field.push_back(rfv);
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::computeRelativeNormalField(const Array<Real> & field,
					    Array<Real> & rel_normal_field) const {
  AKANTU_DEBUG_IN();
  
  // resize arrays to zero
  rel_normal_field.resize(0);
  
  UInt dim = this->model.getSpatialDimension();
  Real * field_p = field.storage();
  Real * normals_p = this->normals.storage();

  Array<Real>::const_iterator< Vector<Real> > it_field  = field.begin(dim);
  Array<Real>::const_iterator< Vector<Real> > it_normal = this->normals.getArray().begin(dim);

  Vector<Real> rfv(dim);

  UInt nb_contact_nodes = this->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt slave  = this->slaves(n);
    UInt master = this->masters(n);

    // relative field vector (slave - master) 
    rfv = it_field[slave];
    rfv -= it_field[master];

    // length of normal projection of relative field
    const Vector<Real> & normal_v = it_normal[n];
    rel_normal_field.push_back(rfv.dot(normal_v));
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Int NTNContact::getNodeIndex(UInt node) const {
  AKANTU_DEBUG_IN();

  Int slave_i  = NTNBaseContact::getNodeIndex(node);
  Int master_i = this->masters.find(node);

  AKANTU_DEBUG_OUT();
  return std::max(slave_i,master_i);
}

/* -------------------------------------------------------------------------- */
void NTNContact::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNContact [" << std::endl;
  NTNBaseContact::printself(stream, indent);
  stream << space << " + masters       : " << std::endl;
  this->masters.printself(stream, indent + 2);
  stream << space << " + lumped_boundary_mastres : " << std::endl;
  this->lumped_boundary_masters.printself(stream, indent + 2);

  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::syncArrays(SyncChoice sync_choice) {
  AKANTU_DEBUG_IN();
  
  NTNBaseContact::syncArrays(sync_choice);
  
  this->masters.syncElements(sync_choice);
  this->lumped_boundary_masters.syncElements(sync_choice);
    
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addDumpFieldToDumper(const std::string & dumper_name,
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
  
  if(field_id == "lumped_boundary_master") {
    internalAddDumpFieldToDumper(dumper_name,
				 field_id,
				 new DumperIOHelper::NodalField<Real>(this->lumped_boundary_masters.getArray()));
  }
  else {
    NTNBaseContact::addDumpFieldToDumper(dumper_name, field_id);
  }
  
  /*
#undef ADD_FIELD
#endif
  */

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
