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
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_EXPONENTIAL_HH__
#define __AKANTU_PHASEFIELD_EXPONENTIAL_HH__

namespace akantu {
class PhaseFieldExponential : public PhaseField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PhaseFieldExponential(PhaseFieldModel & model, const ID & id = "");

  ~PhaseFieldExponential() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  void computePhiOnQuad(const Matrix<Real> & /*strain_quad*/,
                        Real & /*phi_quad*/, Real & /*phi_hist_quad*/);

  void computeDrivingForce(ElementType /*el_type*/,
                           GhostType /*ghost_type*/) override;

  inline void computeDrivingForceOnQuad(const Real & /*phi_quad*/,
                                        Real & /*driving_force_quad*/);

  inline void computeDamageEnergyDensityOnQuad(const Real & /*phi_quad*/,
                                               Real & /*dam_energy_quad*/,
                                               const Real & /*g_c_quad*/);

public:
  void updateInternalParameters() override;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#include "phasefield_exponential_inline_impl.hh"

#endif
