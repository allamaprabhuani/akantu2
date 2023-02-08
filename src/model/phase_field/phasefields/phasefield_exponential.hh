/**
 * @file   phasefield_exponential.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 19 2020
 * @date last modification: Wed Jun 23 2021
 *
 * @brief  Phasefield law for approximating discrete crack as an exponential
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
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
public:
  /// compute the dissiapted energy
  void computeDissipatedEnergy(ElementType el_type) override;

  void
  computeDissipatedEnergyByElement(ElementType type, UInt index,
                                   Vector<Real> & epot_on_quad_points) override;

protected:
  void computePhiOnQuad(const Matrix<Real> & /*strain_quad*/,
                        Real & /*phi_quad*/, Real & /*phi_hist_quad*/);

  void computeDrivingForce(const ElementType & /*el_type*/,
                           GhostType /*ghost_type*/) override;

  inline void computeDrivingForceOnQuad(const Real & /*phi_quad*/,
                                        Real & /*driving_force_quad*/);

  inline void computeDamageEnergyDensityOnQuad(const Real & /*phi_quad*/,
                                               Real & /*dam_energy_quad*/,
                                               const Real & /*g_c_quad*/);

  inline void
  computeDissipatedEnergyOnQuad(const Real & /*dam_quad*/,
                                const Vector<Real> & /*grad_d_quad */,
                                Real & /*energy*/, Real & /*g_c_quad*/);

public:
  void updateInternalParameters() override;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#include "phasefield_exponential_inline_impl.hh"

#endif
