/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_plastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_VISCOPLASTIC_HH_
#define AKANTU_MATERIAL_VISCOPLASTIC_HH_

namespace akantu {

/**
 * Material plastic isotropic
 *
 * parameters in the material files :
 *   - h : Hardening parameter (default: 0)
 *   - sigmay : Yield stress
 *   - rate : Rate sensitivity
 *   - edot0 : Reference strain rate
 *
 *   - ts: Time step
 */

template <Int spatial_dimension>
class MaterialViscoPlastic : public MaterialPlastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialViscoPlastic(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

protected:
  inline void
  computeStressOnQuad(const Matrix<Real> & grad_u,
                      const Matrix<Real> & previous_grad_u,
                      Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
                      Matrix<Real> & inelastic_strain,
                      const Matrix<Real> & previous_inelastic_strain,
                      Real & iso_hardening) const;

  inline void computeTangentModuliOnQuad(
      Matrix<Real> & tangent, const Matrix<Real> & grad_u,
      const Matrix<Real> & previous_grad_u, const Matrix<Real> & sigma_tensor,
      const Matrix<Real> & previous_sigma_tensor,
      const Real & iso_hardening) const;
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Rate sensitivity component (rate)
  Real rate;

  /// Reference strain rate (edot0)
  Real edot0;

  /// Time step (ts)
  Real ts;
};
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_viscoplastic_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_VISCOPLASTIC_HH_ */
