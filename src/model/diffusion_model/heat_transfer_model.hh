/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "diffusion_model.hh"

#ifndef AKANTU_HEAT_TRANSFER_MODEL_HH_
#define AKANTU_HEAT_TRANSFER_MODEL_HH_

namespace akantu {

class HeatTransferModel : public DiffusionModel {
public:
  HeatTransferModel(Mesh & mesh, Int dim = _all_dimensions,
                    const ID & id = "heat_transfer",
                    const std::shared_ptr<DOFManager> & dof_manager = nullptr)
      : DiffusionModel(mesh, dim, id, dof_manager, "temperature",
                       ModelType::_heat_transfer_model) {}

  [[nodiscard]] Array<Real> & getTemperature() { return this->getDiffusion(); }
  [[nodiscard]] const Array<Real> & getTemperature() const {
    return this->getDiffusion();
  }
  [[nodiscard]] Int getTemperatureRelease() const {
    return this->getDiffusionRelease();
  }
  [[nodiscard]] const Array<Real> & getTemperatureRate() const {
    return this->getDiffusionRate();
  }
  [[nodiscard]] Array<Real> & getExternalHeatRate() {
    return this->getExternalFlow();
  }
  [[nodiscard]] const Array<Real> & getExternalHeatRate() const {
    return this->getExternalFlow();
  }
  [[nodiscard]] const Array<Real> & getInternalHeatRate() const {
    return this->getInternalFlow();
  }
};

} // namespace akantu

#endif /* AKANTU_HEAT_TRANSFER_MODEL_HH_ */
