/**
 * Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "phasefield.hh"
#include "aka_common.hh"
#include "phase_field_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseField::PhaseField(PhaseFieldModel & model, const ID & id)
    : Parsable(ParserType::_phasefield, id), id(id), fem(model.getFEEngine()),
      model(model), g_c("g_c", *this),
      spatial_dimension(this->model.getSpatialDimension()),
      element_filter("element_filter", id), damage_on_qpoints("damage", *this),
      gradd("grad_d", *this), phi("phi", *this), strain("strain", *this),
      driving_force("driving_force", *this),
      driving_energy("driving_energy", *this),
      damage_energy("damage_energy", *this),
      damage_energy_density("damage_energy_density", *this),
      dissipated_energy("dissipated_energy", *this) {

  AKANTU_DEBUG_IN();

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(model.getMesh(),
                            _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PhaseField::PhaseField(PhaseFieldModel & model, Int dim, const Mesh & mesh,
                       FEEngine & fe_engine, const ID & id)
    : Parsable(ParserType::_phasefield, id), id(id), fem(fe_engine),
      model(model), g_c("g_c", *this),
      spatial_dimension(this->model.getSpatialDimension()),
      element_filter("element_filter", id),
      damage_on_qpoints("damage", *this, dim, fe_engine, this->element_filter),
      gradd("grad_d", *this, dim, fe_engine, this->element_filter),
      phi("phi", *this, dim, fe_engine, this->element_filter),
      strain("strain", *this, dim, fe_engine, this->element_filter),
      driving_force("driving_force", *this, dim, fe_engine,
                    this->element_filter),
      driving_energy("driving_energy", *this, dim, fe_engine,
                     this->element_filter),
      damage_energy("damage_energy", *this, dim, fe_engine,
                    this->element_filter),
      damage_energy_density("damage_energy_density", *this, dim, fe_engine,
                            this->element_filter),
      dissipated_energy("dissipated_energy", *this, dim, fe_engine,
                        this->element_filter) {

  AKANTU_DEBUG_IN();

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(mesh, _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PhaseField::~PhaseField() = default;

/* -------------------------------------------------------------------------- */
void PhaseField::initialize() {
  registerParam("l0", l0, Real(0.), _pat_parsable | _pat_readable,
                "length scale parameter");
  registerParam("gc", g_c, _pat_parsable | _pat_readable,
                "critical local fracture energy density");
  registerParam("E", E, _pat_parsable | _pat_readable, "Young's modulus");
  registerParam("nu", nu, _pat_parsable | _pat_readable, "Poisson ratio");

  damage_on_qpoints.initialize(1);

  phi.initialize(1);
  driving_force.initialize(1);

  driving_energy.initialize(spatial_dimension);
  gradd.initialize(spatial_dimension);
  g_c.initialize(1);

  strain.initialize(spatial_dimension * spatial_dimension);

  dissipated_energy.initialize(1);

  damage_energy_density.initialize(1);
  damage_energy.initialize(spatial_dimension * spatial_dimension);
}

/* -------------------------------------------------------------------------- */
void PhaseField::initPhaseField() {
  AKANTU_DEBUG_IN();

  this->phi.initializeHistory();

  this->resizeInternals();

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::resizeInternals() {
  AKANTU_DEBUG_IN();
  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    it->second->resize();
  }

  for (auto it = internal_vectors_int.begin(); it != internal_vectors_int.end();
       ++it) {
    it->second->resize();
  }

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it) {
    it->second->resize();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::updateInternalParameters() {
  this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - 2 * this->nu));
  this->mu = this->E / (2 * (1 + this->nu));
}

/* -------------------------------------------------------------------------- */
void PhaseField::computeAllDrivingForces(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  Int spatial_dimension = model.getSpatialDimension();
  auto & damage = model.getDamage();

  for (const auto & type :
       element_filter.elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_filter = element_filter(type, ghost_type);

    if (elem_filter.empty()) {
      continue;
    }

    // compute the damage on quadrature points
    auto & damage_interpolated = damage_on_qpoints(type, ghost_type);
    fem.interpolateOnIntegrationPoints(damage, damage_interpolated, 1, type,
                                       ghost_type);

    auto & gradd_vect = gradd(type, _not_ghost);
    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(damage, gradd_vect, 1, type, ghost_type,
                                    elem_filter);

    computeDrivingForce(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::assembleInternalForces(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  Array<Real> & internal_force = model.getInternalForce();

  for (auto type : element_filter.elementTypes(_ghost_type = ghost_type)) {
    auto & elem_filter = element_filter(type, ghost_type);
    if (elem_filter.empty()) {
      continue;
    }

    auto nb_element = elem_filter.size();
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    auto & driving_force_vect = driving_force(type, ghost_type);

    Array<Real> nt_driving_force(nb_quadrature_points, nb_nodes_per_element);
    fem.computeNtb(driving_force_vect, nt_driving_force, type, ghost_type,
                   elem_filter);

    Array<Real> int_nt_driving_force(nb_element * nb_quadrature_points,
                                     nb_nodes_per_element);
    fem.integrate(nt_driving_force, int_nt_driving_force, nb_nodes_per_element,
                  type, ghost_type, elem_filter);

    // damage_energy_on_qpoints = gc*l0 = scalar
    auto & driving_energy_vect = driving_energy(type, ghost_type);

    Array<Real> bt_driving_energy(nb_element * nb_quadrature_points,
                                  nb_nodes_per_element);
    fem.computeBtD(driving_energy_vect, bt_driving_energy, type, ghost_type,
                   elem_filter);

    Array<Real> int_bt_driving_energy(nb_element, nb_nodes_per_element);
    fem.integrate(bt_driving_energy, int_bt_driving_energy,
                  nb_nodes_per_element, type, ghost_type, elem_filter);

    model.getDOFManager().assembleElementalArrayLocalArray(
        int_nt_driving_force, internal_force, type, ghost_type, -1,
        elem_filter);

    model.getDOFManager().assembleElementalArrayLocalArray(
        int_bt_driving_energy, internal_force, type, ghost_type, -1,
        elem_filter);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_filter = element_filter(type, ghost_type);
    if (elem_filter.empty()) {
      AKANTU_DEBUG_OUT();
      return;
    }

    auto nb_element = elem_filter.size();
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    auto nt_b_n = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "N^t*b*N");

    auto bt_d_b = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    // damage_energy_density_on_qpoints = gc/l0 + phi = scalar
    auto & damage_energy_density_vect = damage_energy_density(type, ghost_type);

    // damage_energy_on_qpoints = gc*l0 = scalar
    auto & damage_energy_vect = damage_energy(type, ghost_type);

    fem.computeBtDB(damage_energy_vect, *bt_d_b, 2, type, ghost_type,
                    elem_filter);

    fem.computeNtbN(damage_energy_density_vect, *nt_b_n, type, ghost_type,
                    elem_filter);

    /// compute @f$ K_{\grad d} = \int_e \mathbf{N}^t * \mathbf{w} *
    /// \mathbf{N}@f$
    auto K_n = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_n");

    fem.integrate(*nt_b_n, *K_n, nb_nodes_per_element * nb_nodes_per_element,
                  type, ghost_type, elem_filter);

    model.getDOFManager().assembleElementalMatricesToMatrix(
        "K", "damage", *K_n, type, _not_ghost, _symmetric, elem_filter);

    /// compute @f$ K_{\grad d} = \int_e \mathbf{B}^t * \mathbf{W} *
    /// \mathbf{B}@f$
    auto K_b = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_b");

    fem.integrate(*bt_d_b, *K_b, nb_nodes_per_element * nb_nodes_per_element,
                  type, ghost_type, elem_filter);

    model.getDOFManager().assembleElementalMatricesToMatrix(
        "K", "damage", *K_b, type, _not_ghost, _symmetric, elem_filter);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::computeDissipatedEnergyByElements() {
  AKANTU_DEBUG_IN();

  const Array<Real> & damage = model.getDamage();

  for (auto type : element_filter.elementTypes(spatial_dimension, _not_ghost)) {

    Array<Idx> & elem_filter = element_filter(type, _not_ghost);
    if (elem_filter.empty()) {
      continue;
    }

    Array<Real> & damage_interpolated = damage_on_qpoints(type, _not_ghost);

    // compute the damage on quadrature points
    fem.interpolateOnIntegrationPoints(damage, damage_interpolated, 1, type,
                                       _not_ghost);

    Array<Real> & gradd_vect = gradd(type, _not_ghost);

    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(damage, gradd_vect, 1, type, _not_ghost,
                                    elem_filter);

    computeDissipatedEnergy(type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::computeDissipatedEnergy(ElementType /*unused*/) {
  AKANTU_DEBUG_IN();
  AKANTU_TO_IMPLEMENT();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real PhaseField::getEnergy() {
  AKANTU_DEBUG_IN();
  Real edis = 0.;

  computeDissipatedEnergyByElements();

  /// integrate the dissipated energy for each type of elements
  for (auto type : element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    edis += fem.integrate(dissipated_energy(type, _not_ghost), type, _not_ghost,
                          element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return edis;
}

/* -------------------------------------------------------------------------- */
Real PhaseField::getEnergy(ElementType type, Idx index) {
  Real edis = 0.;

  Vector<Real> edis_on_quad_points(fem.getNbIntegrationPoints(type));

  computeDissipatedEnergyByElement(type, index, edis_on_quad_points);

  edis = fem.integrate(edis_on_quad_points, Element{type, index, _not_ghost});

  return edis;
}

/* -------------------------------------------------------------------------- */
Real PhaseField::getEnergy(const Element & element) {
  return getEnergy(element.type, element.element);
}

/* -------------------------------------------------------------------------- */
void PhaseField::beforeSolveStep() {
  this->savePreviousState();
  this->computeAllDrivingForces(_not_ghost);
}

/* -------------------------------------------------------------------------- */
void PhaseField::afterSolveStep() {}

/* -------------------------------------------------------------------------- */
void PhaseField::savePreviousState() {
  AKANTU_DEBUG_IN();

  for (auto pair : internal_vectors_real) {
    if (pair.second->hasHistory()) {
      pair.second->saveCurrentValues();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);
  std::string type = getID().substr(getID().find_last_of(':') + 1);

  stream << space << "PhaseField Material " << type << " [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;
}

} // namespace akantu
