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
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_DAMAGE_LINEAR_HH_
#define AKANTU_MATERIAL_DAMAGE_LINEAR_HH_

namespace akantu {

/**
 * Material liner damage
 *
 * parameters in the material files :
 *   - Sigc : (default: 1e5)
 *   - Gc  : (default: 2)
 */
template <Int spatial_dimension>
class MaterialDamageLinear : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDamageLinear(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialDamageLinear(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(Matrix<Real> & F, Matrix<Real> & sigma,
                                  Real & damage, Real & K);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// kind of toughness
  Real Gc;

  /// critical stress
  Real Sigc;

  /// damage internal variable
  InternalField<Real> K;

  Real Epsmin, Epsmax;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_damage_linear_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_DAMAGE_LINEAR_HH_ */
