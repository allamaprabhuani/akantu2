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
protected:
  void computePhiOnQuad(const Matrix<Real> & /*strain_quad*/,
                        Real & /*phi_quad*/, Real & /*phi_hist_quad*/);

  void computeDrivingForce(ElementType /*el_type*/,
                           GhostType /*ghost_type*/) override;

  inline void computeDrivingForceOnQuad(const Real & /*phi_quad*/,
                                        Real & /*driving_force_quad*/);

  inline void computeDamageEnergyDensityOnQuad(const Real & /*phi_quad*/,
                                               Real & /*dam_energy_quad*/);

public:
  void updateInternalParameters() override;
};

/* -------------------------------------------------------------------------- */
inline void
PhaseFieldExponential::computeDrivingForceOnQuad(const Real & phi_quad,
                                                 Real & driving_force_quad) {
  driving_force_quad = 2.0 * phi_quad;
}

/* -------------------------------------------------------------------------- */
inline void PhaseFieldExponential::computeDamageEnergyDensityOnQuad(
    const Real & phi_quad, Real & dam_energy_quad) {
  dam_energy_quad = 2.0 * phi_quad + this->g_c / this->l0;
}

/* -------------------------------------------------------------------------- */
inline void
PhaseFieldExponential::computePhiOnQuad(const Matrix<Real> & strain_quad,
                                        Real & phi_quad, Real & phi_hist_quad) {

  Matrix<Real> strain_plus(spatial_dimension, spatial_dimension);
  // Matrix<Real> strain_minus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_dir(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_plus(spatial_dimension, spatial_dimension);
  // Matrix<Real> strain_diag_minus(spatial_dimension, spatial_dimension);

  Vector<Real> strain_values(spatial_dimension);

  Real trace_plus;

  strain_plus.zero();
  // strain_minus.zero();
  strain_dir.zero();
  strain_values.zero();
  strain_diag_plus.zero();
  // strain_diag_minus.zero();

  strain_quad.eigh(strain_values, strain_dir);

  for (Int i = 0; i < spatial_dimension; i++) {
    strain_diag_plus(i, i) = std::max(Real(0.), strain_values(i));
    // strain_diag_minus(i, i) = std::min(Real(0.), strain_values(i));
  }

  strain_plus = strain_dir * strain_diag_plus * strain_dir.transpose();

  // Nicolas: @mohit is the second transpose really here ?
  // strain_minus = strain_dir * strain_diag_minus * strain_dir.transpose();

  trace_plus = std::max(Real(0.), strain_quad.trace());
  //  trace_minus = std::min(Real(0.), strain_quad.trace());

  auto I = Matrix<Real>::Identity(spatial_dimension, spatial_dimension);

  Matrix<Real> sigma_plus = I * lambda * trace_plus + 2 * mu * strain_plus;
  // Matrix<Real> sigma_minus = I * lambda * trace_minus + 2 * mu *
  // strain_minus;

  phi_quad = 0.5 * sigma_plus.doubleDot(strain_quad);
  if (phi_quad < phi_hist_quad) {
    phi_quad = phi_hist_quad;
  }
}

} // namespace akantu

#endif
