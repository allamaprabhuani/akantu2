/**
 * @file   solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 14:35:38 2010
 *
 * @brief  Implementation of the SolidMechanicsModel class
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
#include "solid_mechanics_model.hh"
#include "aka_math.hh"
#include "integration_scheme_2nd_order.hh"

#include "static_communicator.hh"
#include "sparse_matrix.hh"
#include "solver.hh"

#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::SolidMechanicsModel(Mesh & mesh,
					 UInt dim,
					 const ModelID & id,
					 const MemoryID & memory_id) :
  Model(id, memory_id),
  time_step(NAN), f_m2a(1.0),
  stiffness_matrix(NULL),
  integrator(new CentralDifference()),
  increment_flag(false), solver(NULL),
  spatial_dimension(dim) {
  AKANTU_DEBUG_IN();

  if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();
  registerFEMObject<MyFEMType>("SolidMechanicsFEM",mesh,spatial_dimension);

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
  registerTag(_gst_smm_for_strain, "Explicit Residual");
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

  if(solver) delete solver;
  if(stiffness_matrix) delete stiffness_matrix;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Initialisation                                                             */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initVectors() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEM().getMesh().getNbNodes();
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

  const Mesh::ConnectivityTypeList & type_list = getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    UInt nb_element           = getFEM().getMesh().getNbElement(*it);

    if(!element_material[*it]) {
      std::stringstream sstr_elma; sstr_elma << id << ":element_material:" << *it;
      element_material[*it] = &(alloc<UInt>(sstr_elma.str(), nb_element, 1, 0));
    }
  }

  const Mesh::ConnectivityTypeList & ghost_type_list =
    getFEM().getMesh().getConnectivityTypeList(_ghost);

  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    UInt nb_element           = getFEM().getMesh().getNbGhostElement(*it);

    std::stringstream sstr_elma; sstr_elma << id << ":ghost_element_material:" << *it;
    ghost_element_material[*it] = &(alloc<UInt>(sstr_elma.str(), nb_element, 1, 0));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initModel() {
  /// \todo add  the current position  as a parameter to  initShapeFunctions for
  /// large deformation
  getFEM().initShapeFunctions(_not_ghost);

  getFEM().initShapeFunctions(_ghost);
}



/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = getFEM().getMesh().getNbNodes();

  current_position->resize(nb_nodes);
  //Vector<Real> * current_position = new Vector<Real>(nb_nodes, spatial_dimension, NAN, "position");
  Real * current_position_val = current_position->values;
  Real * position_val         = getFEM().getMesh().getNodes().values;
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
  UInt nb_nodes = getFEM().getMesh().getNbNodes();
  residual->resize(nb_nodes);

  /// copy the forces in residual for boundary conditions
  memcpy(residual->values, force->values, nb_nodes*spatial_dimension*sizeof(Real));

  updateCurrentPosition();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* Explicit scheme                                                            */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateResidual(bool need_initialize) {
  AKANTU_DEBUG_IN();

  if (need_initialize) initializeUpdateResidualData();

  /// start synchronization
  asynchronousSynchronize(_gst_smm_for_strain);

  /// call update residual on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->updateResidual(*current_position, _not_ghost);
    //    (*mat_it)->updateResidual(*displacement, _not_ghost);
  }

  /// finalize communications
  waitEndSynchronize(_gst_smm_for_strain);

  /// call update residual on each ghost elements
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->updateResidual(*current_position, _ghost);
    //    (*mat_it)->updateResidual(*displacement, _ghost);
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
    memcpy(increment->values,
	   displacement->values,
	   displacement->getSize()*displacement->getNbComponent()*sizeof(Real));
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
/* Implicit scheme                                                            */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initImplicitSolver() {
  AKANTU_DEBUG_IN();

  std::stringstream sstr; sstr << id << ":stiffness_matrix";

  stiffness_matrix = new SparseMatrix(getFEM().getMesh(), _symmetric,
				      spatial_dimension, sstr.str(), memory_id);

  stiffness_matrix->buildProfile();

#ifdef AKANTU_USE_MUMPS
  // std::stringstream sstr_solv; sstr_solv << id << ":solver_stiffness_matrix";
  // solver = new SolverMumps(*stiffness_matrix, sstr_solv.str());

  // solver->initialize();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  updateCurrentPosition();

  /// start synchronization
  asynchronousSynchronize(_gst_smm_for_strain);

  stiffness_matrix->clear();

  /// call compute stiffness matrix on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->assembleStiffnessMatrix(*current_position, _not_ghost);
  }

  // UInt nb_nodes = getFEM().getMesh().getNbNodes();
  // residual->resize(nb_nodes);

  // /// copy the forces in residual for boundary conditions
  // memcpy(residual->values, displacement->values, nb_nodes*spatial_dimension*sizeof(Real));

  // *residual *= *stiffness_matrix;

  // Real * residual_val = residual->values;
  // Real * force_val = force->values;

  // for (UInt n = 0; n < spatial_dimension*nb_nodes; ++n) {
  //   *residual_val = *force_val - *residual_val;
  //   force_val++; residual_val++;
  // }

  // /// finalize communications
  // waitEndSynchronize(_gst_smm_for_strain);

  // /// call compute stiffness matrix on each ghost elements
  // for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
  //   (*mat_it)->computeStiffnessMatrix(*current_position, _ghost);
  // }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::solve() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving an implicit step.");

  //  stiffness_matrix->applyBoundary(*boundary);
  stiffness_matrix->removeBoundary(*boundary);
  stiffness_matrix->saveMatrix("K.mtx");

  UInt nb_nodes = displacement->getSize();
  UInt nb_degre_of_freedom = displacement->getNbComponent();

  Vector<Real> * tmp = new Vector<Real>(0, 1);
  for (UInt i = 0; i < nb_degre_of_freedom * nb_nodes; ++i) {
    if(! boundary->values[i]) {
      tmp->push_back(residual->values[i]);
    }
  }

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solv; sstr_solv << id << ":solver_stiffness_matrix";
  solver = new SolverMumps(*stiffness_matrix, sstr_solv.str());
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  solver->initialize();

  //  solver->setRHS(*residual);
  solver->setRHS(*tmp);

  if(!increment) setIncrementFlagOn();

  tmp->clear();

  //  solver->solve(*increment);
  solver->solve(*tmp);

  Real * increment_val     = increment->values;
  Real * displacement_val  = displacement->values;
  bool * boundary_val      = boundary->values;
  Real * tmp_val           = tmp->values;

  for (UInt n = 0; n < nb_nodes * nb_degre_of_freedom; ++n) {
    if(!(*boundary_val)) {
      *increment_val = *(tmp_val++);
      *displacement_val += *increment_val;
    }

    displacement_val++;
    boundary_val++;
    increment_val++;
  }

  stiffness_matrix->restoreProfile();

  delete tmp;
  delete solver;
  solver = NULL;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceIncrement(Real tolerance) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = displacement->getSize();
  UInt nb_degre_of_freedom = displacement->getNbComponent();

  Real norm = 0;
  Real * increment_val     = increment->values;
  bool * boundary_val      = boundary->values;

  for (UInt n = 0; n < nb_nodes * nb_degre_of_freedom; ++n) {
    if(!(*boundary_val)) {
      norm += *increment_val * *increment_val;
    }
    boundary_val++;
    increment_val++;
  }

  AKANTU_DEBUG_INFO("Norm of increment : " << sqrt(norm));

  AKANTU_DEBUG_ASSERT(!isnan(norm), "Something goes wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  return (sqrt(norm) < tolerance);
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceResidual(Real tolerance) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = residual->getSize();

  Real norm = 0;
  Real * residual_val = residual->values;
  bool * boundary_val = boundary->values;

  for (UInt n = 0; n < nb_nodes * spatial_dimension; ++n) {
    if(!(*boundary_val)) {
      norm += *residual_val * *residual_val;
    }
    boundary_val++;
    residual_val++;
  }

  AKANTU_DEBUG_INFO("Norm of residual : " << sqrt(norm));

  AKANTU_DEBUG_ASSERT(!isnan(norm), "Something goes wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  return (sqrt(norm) < tolerance);
}


/* -------------------------------------------------------------------------- */
/* Information                                                                */
/* -------------------------------------------------------------------------- */

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
    UInt nb_nodes = getFEM().getMesh().getNbNodes();
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

  Real * coord    = getFEM().getMesh().getNodes().values;
  Real * disp_val = displacement->values;

  const Mesh::ConnectivityTypeList & type_list = getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(getFEM().getMesh().getSpatialDimension(*it) != spatial_dimension) continue;


    UInt nb_nodes_per_element = getFEM().getMesh().getNbNodesPerElement(*it);
    UInt nb_element           = getFEM().getMesh().getNbElement(*it);

    UInt * conn         = getFEM().getMesh().getConnectivity(*it).values;
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

      Real el_size    = getFEM().getElementInradius(u, *it);
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

  UInt nb_nodes = getFEM().getMesh().getNbNodes();
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
  getFEM().printself(stream, indent + 2);
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
  const Mesh::ConnectivityTypeList & type_list = getFEM().getMesh().getConnectivityTypeList();
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
