/**
 * @file   solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 14:35:38 2010
 *
 * @brief  Implementation of the SolidMechanicsModel class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "aka_math.hh"
#include "integration_scheme/central_difference.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::SolidMechanicsModel(UInt spatial_dimension,
					 const ModelID & id,
					 const MemoryID & memory_id) :
  Model(spatial_dimension, id, memory_id),
  time_step(NAN), f_m2a(1.0),
  integrator(new CentralDifference()) {
  AKANTU_DEBUG_IN();

  this->displacement = NULL;
  this->mass         = NULL;
  this->velocity     = NULL;
  this->acceleration = NULL;
  this->force        = NULL;
  this->residual     = NULL;
  this->boundary     = NULL;

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->stress          [t] = NULL;
    this->strain          [t] = NULL;
    this->element_material[t] = NULL;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::SolidMechanicsModel(Mesh & mesh,
					 UInt spatial_dimension,
					 const ModelID & id,
					 const MemoryID & memory_id) :
  Model(mesh, spatial_dimension, id, memory_id),
  time_step(NAN), f_m2a(1.0),
  integrator(new CentralDifference()) {
  AKANTU_DEBUG_IN();

  this->displacement = NULL;
  this->mass         = NULL;
  this->velocity     = NULL;
  this->acceleration = NULL;
  this->force        = NULL;
  this->residual     = NULL;
  this->boundary     = NULL;


  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->stress          [t] = NULL;
    this->strain          [t] = NULL;
    this->element_material[t] = NULL;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::~SolidMechanicsModel() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG(dblAccessory, "Deleting displacements vector of type");
  dealloc(displacement->getID());
  displacement = NULL;

  AKANTU_DEBUG(dblAccessory, "Deleting mass vector of type");
  dealloc(mass->getID());
  mass = NULL;

  AKANTU_DEBUG(dblAccessory, "Deleting velocity vector of type");
  dealloc(velocity->getID());
  velocity = NULL;

  AKANTU_DEBUG(dblAccessory, "Deleting acceleration vector of type");
  dealloc(acceleration->getID());
  acceleration = NULL;

  AKANTU_DEBUG(dblAccessory, "Deleting force vector of type");
  dealloc(force->getID());
  force = NULL;

  AKANTU_DEBUG(dblAccessory, "Deleting residual vector of type");
  dealloc(residual->getID());
  residual = NULL;

  AKANTU_DEBUG(dblAccessory, "Deleting boundary vector of type");
  dealloc(boundary->getID());
  boundary = NULL;

  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    AKANTU_DEBUG(dblAccessory, "Deleting stress derivatives vector of type " << *it);
    dealloc(stress[*it]->getID());
    stress[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting strain vector of type " << *it);
    dealloc(strain[*it]->getID());
    strain[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting element_material vector of type " << *it);
    dealloc(element_material[*it]->getID());
    element_material[*it] = NULL;
  }

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    delete *mat_it;
  }
  materials.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initVectors() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = fem->getNbNodes();
  std::stringstream sstr_disp; sstr_disp << id << ":displacement";
  std::stringstream sstr_mass; sstr_mass << id << ":mass";
  std::stringstream sstr_velo; sstr_velo << id << ":velocity";
  std::stringstream sstr_acce; sstr_acce << id << ":acceleration";
  std::stringstream sstr_forc; sstr_forc << id << ":force";
  std::stringstream sstr_resi; sstr_resi << id << ":residual";
  std::stringstream sstr_boun; sstr_boun << id << ":boundary";

#ifdef AKANTU_DEBUG
  Real init_val = NAN;
#else
  Real init_val = 0;
#endif

  displacement = &(alloc<Real>(sstr_disp.str(), nb_nodes, spatial_dimension, init_val));
  mass         = &(alloc<Real>(sstr_mass.str(), nb_nodes, 1)); // \todo see if it must not be spatial_dimension
  velocity     = &(alloc<Real>(sstr_velo.str(), nb_nodes, spatial_dimension, init_val));
  acceleration = &(alloc<Real>(sstr_acce.str(), nb_nodes, spatial_dimension, init_val));
  force        = &(alloc<Real>(sstr_forc.str(), nb_nodes, spatial_dimension, init_val));
  residual     = &(alloc<Real>(sstr_resi.str(), nb_nodes, spatial_dimension, init_val));
  boundary     = &(alloc<bool> (sstr_boun.str(), nb_nodes, spatial_dimension, false));

  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    UInt nb_quadrature_points = fem->getNbQuadraturePoints(*it);
    UInt nb_element           = fem->getNbElement(*it);

    std::stringstream sstr_stre; sstr_stre << id << ":stress:" << *it;
    std::stringstream sstr_stra; sstr_stra << id << ":strain:" << *it;
    std::stringstream sstr_elma; sstr_elma << id << ":element_material:" << *it;
    stress          [*it] = &(alloc<Real>(sstr_stre.str(), nb_element,
					  nb_quadrature_points * spatial_dimension * spatial_dimension,
					  init_val));
    strain          [*it] = &(alloc<Real>(sstr_stra.str(), nb_element,
					  nb_quadrature_points * spatial_dimension * spatial_dimension,
					  init_val));
    element_material[*it] = &(alloc<UInt>(sstr_elma.str(), nb_element, 1, 0));
  }

  AKANTU_DEBUG_OUT();
};

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::readMaterials(const std::string & filename) {
  std::stringstream sstr; sstr << id << ":0:elastic";
  materials.push_back(new MaterialElastic(*this, sstr.str()));
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initMaterials() {
  Material ** mat_val = &(materials.at(0));

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    /// init internals properties
    (*mat_it)->initMaterial();
  }

  /// fill the element filters of the materials using the element_material arrays
  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(fem->getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_element = fem->getNbElement(*it);
    UInt * elem_mat_val = element_material[*it]->values;

    for (UInt el = 0; el < nb_element; ++el) {
      mat_val[elem_mat_val[el]]->element_filter[*it]->push_back(el);
    }
  }
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initModel() {
  fem->initShapeFunctions();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass() {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));

  UInt nb_nodes = fem->getNbNodes();
  memset(mass->values, 0, nb_nodes*sizeof(Real));

  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(fem->getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_nodes_per_element = fem->getNbNodesPerElement(*it);
    UInt nb_element           = fem->getNbElement(*it);

    const Vector<Real> & shapes = fem->getShapes(*it);
    Vector<Real> * rho_phi_i = new Vector<Real>(nb_element, nb_nodes_per_element, "rho_x_shapes");

    UInt * elem_mat_val = element_material[*it]->values;
    Real * rho_phi_i_val = rho_phi_i->values;
    Real * shapes_val = shapes.values;

    /// compute rho * \phi_i for each nodes of each element
    for (UInt el = 0; el < nb_element; ++el) {
      Real rho = mat_val[elem_mat_val[el]]->getRho();
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	*rho_phi_i_val++ = rho * *shapes_val++;
      }
    }

    Vector<Real> * int_rho_phi_i = new Vector<Real>(nb_element, nb_nodes_per_element, "inte_rho_x_shapes");
    fem->integrate(*rho_phi_i, *int_rho_phi_i, nb_nodes_per_element, *it);
    delete rho_phi_i;

    fem->assembleVector(*int_rho_phi_i, *mass, 1, *it);
    delete int_rho_phi_i;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateResidual() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = fem->getNbNodes();

  Vector<Real> * current_position = new Vector<Real>(nb_nodes, spatial_dimension, NAN, "position");
  Real * current_position_val = current_position->values;
  Real * position_val         = fem->getMesh().getNodes().values;
  Real * displacement_val     = displacement->values;

  /// compute current_position = initial_position + displacement
  memcpy(current_position_val, position_val, nb_nodes*spatial_dimension*sizeof(Real));
  for (UInt n = 0; n < nb_nodes*spatial_dimension; ++n) {
    *current_position_val++ += *displacement_val++;
  }

  /// copy the forces in residual for boundary conditions
  memcpy(residual->values, force->values, nb_nodes*spatial_dimension*sizeof(Real));

  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(fem->getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_nodes_per_element = fem->getNbNodesPerElement(*it);
    UInt nb_quadrature_points = fem->getNbQuadraturePoints(*it);
    UInt nb_element           = fem->getNbElement(*it);

    fem->gradientOnQuadraturePoints(*current_position, *strain[*it], spatial_dimension, *it);

    /// compute the constitutive law for each materials
    std::vector<Material *>::iterator mat_it;
    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      (*mat_it)->constitutiveLaw(*it);
    }

    Vector<Real> * sigma_dphi_dx =
      new Vector<Real>(nb_element, nb_nodes_per_element * spatial_dimension * nb_quadrature_points);

    const Vector<Real> & shapes_derivatives = fem->getShapesDerivatives(*it);

    Real * shapesd_val       = shapes_derivatives.values;
    Real * stress_val        = stress[*it]->values;
    Real * sigma_dphi_dx_val = sigma_dphi_dx->values;

    UInt offset_shapesd_val       = spatial_dimension * nb_nodes_per_element;
    UInt offset_stress_val        = spatial_dimension * spatial_dimension;
    UInt offset_sigma_dphi_dx_val = spatial_dimension * nb_nodes_per_element;

    /// compute \sigma * \partial \phi / \partial X
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quadrature_points; ++q) {
	Math::matrix_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
			     shapesd_val, stress_val, sigma_dphi_dx_val);
	shapesd_val       += offset_shapesd_val;
	stress_val        += offset_stress_val;
	sigma_dphi_dx_val += offset_sigma_dphi_dx_val;
      }
    }

    Vector<Real> * int_sigma_dphi_dx = new Vector<Real>(nb_element, nb_nodes_per_element*spatial_dimension);
    fem->integrate(*sigma_dphi_dx, *int_sigma_dphi_dx, nb_nodes_per_element*spatial_dimension, *it);
    delete sigma_dphi_dx;

    fem->assembleVector(*int_sigma_dphi_dx, *residual, residual->getNbComponent(), *it, NULL, -1);
    delete int_sigma_dphi_dx;
  }
  delete current_position;

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

  integrator->integrationSchemePred(time_step,
				    *displacement,
				    *velocity,
				    *acceleration,
				    *boundary);

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
Real SolidMechanicsModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));
  Real min_dt = HUGE_VAL;

  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(fem->getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_element = fem->getNbElement(*it);
    UInt * elem_mat_val = element_material[*it]->values;
    for (UInt el = 0; el < nb_element; ++el) {
      Real el_size    = fem->getElementInradius(el, *it);
      Real el_dt      = mat_val[elem_mat_val[el]]->getStableTimeStep(el_size);
      min_dt = min_dt > el_dt ? el_dt : min_dt;
    }
  }

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setPotentialEnergyFlagOn() {
  AKANTU_DEBUG_IN();
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->setPotentialEnergyFlagOn();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setPotentialEnergyFlagOff() {
  AKANTU_DEBUG_IN();
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->setPotentialEnergyFlagOff();
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  /// fill the element filters of the materials using the element_material arrays
  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(fem->getSpatialDimension(*it) != spatial_dimension) continue;

    std::vector<Material *>::iterator mat_it;
    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      epot += fem->integrate((*mat_it)->getPotentialEnergy(*it),
			     *it);
    }
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = velocity->getSize();
  UInt nb_degre_of_freedom = velocity->getNbComponent();

  Real * mass_val = mass->values;
  Real * vel_val  = velocity->values;

  Real ekin = 0.;

  for (UInt n = 0; n < nb_nodes; ++n) {
    Real norm_vel = 0.;
    for (UInt d = 0; d < nb_degre_of_freedom; d++) {
      norm_vel += *vel_val * *vel_val;
      vel_val++;
    }
    ekin += *mass_val * norm_vel;

    mass_val++;
  }

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
  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    stream << space << AKANTU_INDENT << AKANTU_INDENT << " + " << *it <<" [" << std::endl;
    stress          [*it]->printself(stream, indent + 3);
    strain          [*it]->printself(stream, indent + 3);
    element_material[*it]->printself(stream, indent + 3);
    stream << space << AKANTU_INDENT << AKANTU_INDENT << "]" << std::endl;
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
