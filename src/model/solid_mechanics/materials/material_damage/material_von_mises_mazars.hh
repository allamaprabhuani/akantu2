/**
 * @file   material_von_mises_mazars.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Dec 04 2020
 *
 * @brief  Mazars damage with Von Misses criteria
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_damage.hh"
#include "material_linear_isotropic_hardening.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VONMISES_MAZARS_HH__
#define __AKANTU_MATERIAL_VONMISES_MAZARS_HH__

namespace akantu {

template <UInt spatial_dimension,
          template <UInt> class Parent = MaterialLinearIsotropicHardening>
class MaterialVonMisesMazars
    : public MaterialDamage<spatial_dimension, Parent> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialVonMisesMazars(SolidMechanicsModel & model, const ID & id = "");
  MaterialVonMisesMazars(SolidMechanicsModel & model, UInt dim,
                         const Mesh & mesh, FEEngine & fe_engine,
                         const ID & id = "");
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
                                  Matrix<Real> & sigma, Real & damage,
                                  Real & Ehat);

  inline void computeDamageAndStressOnQuad(const Matrix<Real> & grad_u,
                                           Matrix<Real> & sigma, Real & damage,
                                           Real & Ehat);

  inline void computeDamageOnQuad(const Real & epsilon_equ,
                                  const Matrix<Real> & sigma,
                                  const Vector<Real> & epsilon_princ,
                                  Real & dam);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// damage threshold
  RandomInternalField<Real> K0;
  /// parameter damage traction 1
  Real At;
  /// parameter damage traction 2
  Real Bt;
  /// parameter damage compression 1
  Real Ac;
  /// parameter damage compression 2
  Real Bc;
  /// parameter for shear
  Real beta;

  /// specify the variable to average false = ehat, true = damage (only valid
  /// for non local version)
  bool damage_in_compute_stress;
};

} // namespace akantu

#include "material_von_mises_mazars_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_VONMISES_MAZARS_HH__ */
