/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "coupler_solid_contact.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template <>
CouplerSolidContactTemplate<SolidMechanicsModel>::CouplerSolidContactTemplate(
    Mesh & mesh, Int dim, const ID & id,
    const std::shared_ptr<DOFManager> & dof_manager)
    : Model(mesh, ModelType::_coupler_solid_contact, dim, id) {

  this->initDOFManager(dof_manager);

  this->mesh.registerDumper<DumperParaview>("coupler_solid_contact", id, true);
  this->mesh.addDumpMeshToDumper("coupler_solid_contact", mesh,
                                 Model::spatial_dimension, _not_ghost,
                                 _ek_regular);
  this->registerDataAccessor(*this);

  solid = std::make_unique<SolidMechanicsModel>(mesh, Model::spatial_dimension,
                                                "solid_mechanics_model",
                                                this->dof_manager);
  contact = std::make_unique<ContactMechanicsModel>(
      mesh, Model::spatial_dimension, "contact_mechanics_model",
      this->dof_manager);
}

/* -------------------------------------------------------------------------- */
template <>
void CouplerSolidContactTemplate<SolidMechanicsModel>::initFullImpl(
    const ModelOptions & options) {
  Model::initFullImpl(options);

  solid->initFull(_analysis_method = this->method);
  contact->initFull(_analysis_method = this->method);
}

} // namespace akantu
