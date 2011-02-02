/**
 * @file   solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 14:35:38 2010
 *
 * @brief  Implementation of the SolidMechanicsModel class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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
#include "solid_mechanics_model.hh"
#include "aka_math.hh"
#include "integration_scheme_2nd_order.hh"

#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::SolidMechanicsModel(Mesh & mesh,
					 UInt spatial_dimension,
					 const ModelID & id,
					 const MemoryID & memory_id) :
  Model(mesh, spatial_dimension, id, memory_id),
  time_step(NAN), f_m2a(1.0),
  integrator(new CentralDifference()),
  increment_flag(false) {
  AKANTU_DEBUG_IN();

  this->displacement = NULL;
  this->mass         = NULL;
  this->velocity     = NULL;
  this->acceleration = NULL;
  this->force        = NULL;
  this->residual     = NULL;
  this->boundary     = NULL;

  this->increment    = NULL;


  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->element_material[t] = NULL;
    this->ghost_element_material[t] = NULL;
  }

  registerTag(_gst_smm_mass, "Mass");
  registerTag(_gst_smm_residual, "Explicit Residual");
  registerTag(_gst_smm_boundary, "Boundary conditions");

  materials.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::~SolidMechanicsModel() {
  AKANTU_DEBUG_IN();

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    delete *mat_it;
  }
  materials.clear();

  delete integrator;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initVectors() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = fem->getMesh().getNbNodes();
  std::stringstream sstr_disp; sstr_disp << id << ":displacement";
  std::stringstream sstr_mass; sstr_mass << id << ":mass";
  std::stringstream sstr_velo; sstr_velo << id << ":velocity";
  std::stringstream sstr_acce; sstr_acce << id << ":acceleration";
  std::stringstream sstr_forc; sstr_forc << id << ":force";
  std::stringstream sstr_resi; sstr_resi << id << ":residual";
  std::stringstream sstr_boun; sstr_boun << id << ":boundary";

  displacement = &(alloc<Real>(sstr_disp.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  mass         = &(alloc<Real>(sstr_mass.str(), nb_nodes, 1)); // \todo see if it must not be spatial_dimension
  velocity     = &(alloc<Real>(sstr_velo.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  acceleration = &(alloc<Real>(sstr_acce.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  force        = &(alloc<Real>(sstr_forc.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  residual     = &(alloc<Real>(sstr_resi.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  boundary     = &(alloc<bool>(sstr_boun.str(), nb_nodes, spatial_dimension, false));

  std::stringstream sstr_curp; sstr_curp << id << ":current_position_tmp";
  current_position = &(alloc<Real>(sstr_curp.str(), 0, spatial_dimension, REAL_INIT_VALUE));

  const Mesh::ConnectivityTypeList & type_list = fem->getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    UInt nb_element           = fem->getMesh().getNbElement(*it);

    if(!element_material[*it]) {
      std::stringstream sstr_elma; sstr_elma << id << ":element_material:" << *it;
      element_material[*it] = &(alloc<UInt>(sstr_elma.str(), nb_element, 1, 0));
    }
  }

  const Mesh::ConnectivityTypeList & ghost_type_list =
    fem->getMesh().getConnectivityTypeList(_ghost);

  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    UInt nb_element           = fem->getMesh().getNbGhostElement(*it);

    std::stringstream sstr_elma; sstr_elma << id << ":ghost_element_material:" << *it;
    ghost_element_material[*it] = &(alloc<UInt>(sstr_elma.str(), nb_element, 1, 0));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initModel() {
  fem->initShapeFunctions(_not_ghost);
  fem->initShapeFunctions(_ghost);
}



/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = fem->getMesh().getNbNodes();

  current_position->resize(nb_nodes);
  //Vector<Real> * current_position = new Vector<Real>(nb_nodes, spatial_dimension, NAN, "position");
  Real * current_position_val = current_position->values;
  Real * position_val         = fem->getMesh().getNodes().values;
  Real * displacement_val     = displacement->values;

  /// compute current_position = initial_position + displacement
  memcpy(current_position_val, position_val, nb_nodes*spatial_dimension*sizeof(Real));
  for (UInt n = 0; n < nb_nodes*spatial_dimension; ++n) {
    *current_position_val++ += *displacement_val++;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initializeUpdateResidualData() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = fem->getMesh().getNbNodes();
  residual->resize(nb_nodes);

  /// copy the forces in residual for boundary conditions
  memcpy(residual->values, force->values, nb_nodes*spatial_dimension*sizeof(Real));

  updateCurrentPosition();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateResidual(bool need_initialize) {
  AKANTU_DEBUG_IN();

  if (need_initialize) initializeUpdateResidualData();

  /// start synchronization
  asynchronousSynchronize(_gst_smm_residual);

  /// call update residual on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->updateResidual(*current_position, _not_ghost);
  }

  /// finalize communications
  waitEndSynchronize(_gst_smm_residual);

  /// call update residual on each ghost elements
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->updateResidual(*current_position, _ghost);
  }

  //  current_position;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateAcceleration() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = acceleration->getSize();

  UInt nb_degre_of_freedom = acceleration->getNbComponent();

  Real * mass_val     = mass->values;
  Real * residual_val = residual->values;
  Real * accel_val    = acceleration->values;
  bool * boundary_val = boundary->values;

  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt d = 0; d < nb_degre_of_freedom; d++) {
      if(!(*boundary_val)) {
	*accel_val = f_m2a * *residual_val / *mass_val;
      }
      residual_val++;
      accel_val++;
      boundary_val++;
    }
    mass_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::explicitPred() {
  AKANTU_DEBUG_IN();

  if(increment_flag) {
    memcpy(increment->values, displacement->values, displacement->getSize()*displacement->getNbComponent()*sizeof(Real));
  }

  integrator->integrationSchemePred(time_step,
				    *displacement,
				    *velocity,
				    *acceleration,
				    *boundary);

  if(increment_flag) {
    Real * inc_val = increment->values;
    Real * dis_val = displacement->values;
    UInt nb_nodes = displacement->getSize();

    for (UInt n = 0; n < nb_nodes; ++n) {
      *inc_val = *dis_val - *inc_val;
      *inc_val++;
      *dis_val++;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::explicitCorr() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemeCorr(time_step,
				    *displacement,
				    *velocity,
				    *acceleration,
				    *boundary);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::synchronizeBoundaries() {
  AKANTU_DEBUG_IN();
  synchronize(_gst_smm_boundary);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setIncrementFlagOn() {
  AKANTU_DEBUG_IN();

  if(!increment) {
    UInt nb_nodes = fem->getMesh().getNbNodes();
    std::stringstream sstr_inc; sstr_inc << id << ":increment";
    increment = &(alloc<Real>(sstr_inc.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  }

  increment_flag = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));
  Real min_dt = std::numeric_limits<Real>::max();

  Real * coord    = fem->getMesh().getNodes().values;
  Real * disp_val = displacement->values;

  const Mesh::ConnectivityTypeList & type_list = fem->getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(fem->getMesh().getSpatialDimension(*it) != spatial_dimension) continue;


    UInt nb_nodes_per_element = fem->getMesh().getNbNodesPerElement(*it);
    UInt nb_element           = fem->getMesh().getNbElement(*it);

    UInt * conn         = fem->getMesh().getConnectivity(*it).values;
    UInt * elem_mat_val = element_material[*it]->values;
    Real * u = new Real[nb_nodes_per_element*spatial_dimension];

    for (UInt el = 0; el < nb_element; ++el) {
      UInt el_offset  = el * nb_nodes_per_element;

      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	UInt offset_conn = conn[el_offset + n] * spatial_dimension;
	memcpy(u + n * spatial_dimension,
	       coord + offset_conn,
	       spatial_dimension * sizeof(Real));

	for (UInt i = 0; i < spatial_dimension; ++i) {
	  u[n * spatial_dimension + i] += disp_val[offset_conn + i];
	}
      }

      Real el_size    = fem->getElementInradius(u, *it);
      Real el_dt      = mat_val[elem_mat_val[el]]->getStableTimeStep(el_size);
      min_dt = min_dt > el_dt ? el_dt : min_dt;
    }

    delete [] u;
  }


  /// reduction min over all processors
  allReduce(&min_dt, _so_min);


  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  /// call update residual on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    epot += (*mat_it)->getPotentialEnergy();
  }

  /// reduction sum over all processors
  allReduce(&epot, _so_sum);

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy() {
  AKANTU_DEBUG_IN();

  Real ekin = 0.;

  UInt nb_nodes = fem->getMesh().getNbNodes();
  //  Vector<Real> * v_square = new Vector<Real>(nb_nodes, 1, "v_square");

  Real * vel_val  = velocity->values;
  Real * mass_val = mass->values;

  for (UInt n = 0; n < nb_nodes; ++n) {
    Real v2 = 0;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      v2 += *vel_val * *vel_val;
      vel_val++;
    }
    ekin += *mass_val * v2;
    mass_val++;
  }


  // Real * v_s_val  = v_square->values;

  // for (UInt n = 0; n < nb_nodes; ++n) {
  //   *v_s_val = 0;
  //   for (UInt s = 0; s < spatial_dimension; ++s) {
  //     *v_s_val += *vel_val * *vel_val;
  //     vel_val++;
  //   }
  //   v_s_val++;
  // }

  // Material ** mat_val = &(materials.at(0));

  // const Mesh:: ConnectivityTypeList & type_list = fem->getMesh().getConnectivityTypeList();
  // Mesh::ConnectivityTypeList::const_iterator it;
  // for(it = type_list.begin(); it != type_list.end(); ++it) {
  //   if(fem->getMesh().getSpatialDimension(*it) != spatial_dimension) continue;

  //   UInt nb_quadrature_points = FEM::getNbQuadraturePoints(*it);
  //   UInt nb_element = fem->getMesh().getNbElement(*it);

  //   Vector<Real> * v_square_el = new Vector<Real>(nb_element * nb_quadrature_points, 1, "v_square per element");

  //   fem->interpolateOnQuadraturePoints(*v_square, *v_square_el, 1, *it);

  //   Real * v_square_el_val = v_square_el->values;
  //   UInt * elem_mat_val = element_material[*it]->values;

  //   for (UInt el = 0; el < nb_element; ++el) {
  //     Real rho = mat_val[elem_mat_val[el]]->getRho();
  //     for (UInt q = 0; q < nb_quadrature_points; ++q) {
  // 	*v_square_el_val *= rho;
  // 	v_square_el_val++;
  //     }
  //   }

  //   ekin += fem->integrate(*v_square_el, *it);

  //   delete v_square_el;
  // }
  // delete v_square;

  /// reduction sum over all processors
  allReduce(&ekin, _so_sum);

  AKANTU_DEBUG_OUT();
  return ekin * .5;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Solid Mechanics Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << spatial_dimension << std::endl;

  stream << space << " + fem [" << std::endl;
  fem->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;
  stream << space << " + nodals information [" << std::endl;
  displacement->printself(stream, indent + 2);
  mass        ->printself(stream, indent + 2);
  velocity    ->printself(stream, indent + 2);
  acceleration->printself(stream, indent + 2);
  force       ->printself(stream, indent + 2);
  residual    ->printself(stream, indent + 2);
  boundary    ->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + connectivity type information [" << std::endl;
  const Mesh::ConnectivityTypeList & type_list = fem->getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    stream << space << AKANTU_INDENT << AKANTU_INDENT << " + " << *it <<" [" << std::endl;
    element_material[*it]->printself(stream, indent + 3);
    stream << space << AKANTU_INDENT << AKANTU_INDENT << "]" << std::endl;
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
