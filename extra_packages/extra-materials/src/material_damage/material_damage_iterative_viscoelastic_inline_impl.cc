/**
 * @file   material_iterative_stiffness_reduction.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Feb 18 16:03:56 2016
 *
 * @brief  Implementation of material iterative stiffness reduction
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
#include "material_damage_iterative_viscoelastic.hh"
#include "communicator.hh"
#include "solid_mechanics_model_RVE.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {


/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, Real & dam) {

  Real dt = this->model.getTimeStep();
  Real E_ef = this->Einf * (1 - dam);

  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Real lambda = this->Eta(k) / this->Ev(k);
    Real exp_dt_lambda = std::exp(-dt / lambda);
    if (exp_dt_lambda == 1) {
      E_ef += (1 - dam) * this->Ev(k);
    } else {
      E_ef += (1 - exp_dt_lambda) * (1 - dam) * this->Ev(k) * lambda / dt;
    }
  }

  tangent.copy(this->C);
  tangent *= E_ef;

}

/* -------------------------------------------------------------------------- */

} // namespace akantu
