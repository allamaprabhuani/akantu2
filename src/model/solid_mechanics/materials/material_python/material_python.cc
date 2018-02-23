/**
 * @file   material_python.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Nov 13 2015
 * @date last modification: Fri Nov 13 2015
 *
 * @brief  Material python implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "material_python.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MaterialPython::MaterialPython(SolidMechanicsModel & model, PyObject * obj,
                               const ID & id)
    : Material(model, id), PythonFunctor(obj) {
  AKANTU_DEBUG_IN();

  this->registerInternals();

  std::vector<std::string> param_names =
      this->callFunctor<std::vector<std::string>>("registerParam");

  for (UInt i = 0; i < param_names.size(); ++i) {
    std::stringstream sstr;
    sstr << "PythonParameter" << i;
    this->registerParam(param_names[i], local_params[param_names[i]], 0.,
                        _pat_parsable | _pat_readable, sstr.str());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialPython::registerInternals() {
  std::vector<std::string> internal_names =
      this->callFunctor<std::vector<std::string>>("registerInternals");

  std::vector<UInt> internal_sizes;

  try {
    internal_sizes =
        this->callFunctor<std::vector<UInt>>("registerInternalSizes");
  } catch (...) {
    internal_sizes.assign(internal_names.size(), 1);
  }

  for (UInt i = 0; i < internal_names.size(); ++i) {
    std::stringstream sstr;
    sstr << "PythonInternal" << i;
    this->internals[internal_names[i]] =
        std::make_unique<InternalField<Real>>(internal_names[i], *this);
    AKANTU_DEBUG_INFO("alloc internal " << internal_names[i] << " "
                                        << &this->internals[internal_names[i]]);

    this->internals[internal_names[i]]->initialize(internal_sizes[i]);
  }

  // making an internal with the quadrature points coordinates
  this->internals["quad_coordinates"] =
      std::make_unique<InternalField<Real>>("quad_coordinates", *this);
  auto && coords = *this->internals["quad_coordinates"];
  coords.initialize(this->getSpatialDimension());
}

/* -------------------------------------------------------------------------- */

void MaterialPython::initMaterial() {
  AKANTU_DEBUG_IN();

  Material::initMaterial();

  auto && coords = *this->internals["quad_coordinates"];
  this->model.getFEEngine().computeIntegrationPointsCoordinates(
      coords, &this->element_filter);

  auto params = local_params;
  params["rho"] = this->rho;

  try {
    this->callFunctor<void>("initMaterial", this->internals, params);
  } catch (...) {
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void MaterialPython::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto params = local_params;
  params["rho"] = this->rho;

  std::map<std::string, Array<Real> *> internal_arrays;
  for (auto & i : this->internals) {
    auto & array = (*i.second)(el_type, ghost_type);
    auto & name = i.first;
    internal_arrays[name] = &array;
  }

  this->callFunctor<void>("computeStress", this->gradu(el_type, ghost_type),
                          this->stress(el_type, ghost_type), internal_arrays,
                          params);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialPython::computeTangentModuli(const ElementType & el_type,
                                          Array<Real> & tangent_matrix,
                                          GhostType ghost_type) {
  auto params = local_params;
  params["rho"] = this->rho;

  std::map<std::string, Array<Real> *> internal_arrays;
  for (auto & i : this->internals) {
    auto & array = (*i.second)(el_type, ghost_type);
    auto & name = i.first;
    internal_arrays[name] = &array;
  }

  this->callFunctor<void>("computeTangentModuli",
                          this->gradu(el_type, ghost_type), tangent_matrix,
                          internal_arrays, params);
}

/* -------------------------------------------------------------------------- */
Real MaterialPython::getPushWaveSpeed(const Element &) const {
  auto params = local_params;
  params["rho"] = this->rho;

  return this->callFunctor<Real>("getPushWaveSpeed", params);
}

/* -------------------------------------------------------------------------- */

Real MaterialPython::getEnergyForType(const std::string & type,
                                      ElementType el_type) {
  AKANTU_DEBUG_IN();

  std::map<std::string, Array<Real> *> internal_arrays;
  for (auto & i : this->internals) {
    auto & array = (*i.second)(el_type, _not_ghost);
    auto & name = i.first;
    internal_arrays[name] = &array;
  }

  auto params = local_params;
  params["rho"] = this->rho;

  auto & energy_density = *internal_arrays[type];

  this->callFunctor<void>("getEnergyDensity", type, energy_density,
                          this->gradu(el_type, _not_ghost),
                          this->stress(el_type, _not_ghost), internal_arrays,
                          params);

  Real energy = fem.integrate(energy_density, el_type, _not_ghost,
                              element_filter(el_type, _not_ghost));

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */

Real MaterialPython::getEnergy(const std::string & type) {
  AKANTU_DEBUG_IN();

  if (this->internals.find(type) == this->internals.end()) {
    AKANTU_EXCEPTION("unknown energy type: "
                     << type << " you must declare an internal named " << type);
  }

  Real energy = 0.;
  /// integrate the potential energy for each type of elements
  Mesh::type_iterator it = element_filter.firstType(spatial_dimension);
  Mesh::type_iterator last_type = element_filter.lastType(spatial_dimension);
  for (; it != last_type; ++it) {
    energy += this->getEnergyForType(type, *it);
  }
  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */

} // akantu
