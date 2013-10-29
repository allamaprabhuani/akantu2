/**
 * @file   material_damage_linear.hh
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Oct 03 11:53:09 2011
 *
 * @brief  Material isotropic elastic + linear softening
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
#include "aka_common.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DAMAGE_LINEAR_HH__
#define __AKANTU_MATERIAL_DAMAGE_LINEAR_HH__

__BEGIN_AKANTU__

/**
 * Material liner damage
 *
 * parameters in the material files :
 *   - Sigc : (default: 1e5)
 *   - Gc  : (default: 2)
 */
template<UInt spatial_dimension>
class MaterialDamageLinear : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialDamageLinear(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialDamageLinear() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(Matrix<Real> & F,
				  Matrix<Real> & sigma,
				  Real & damage,
				  Real &K);

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
  ByElementTypeReal K;

  Real Epsmin, Epsmax;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_damage_linear_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_DAMAGE_LINEAR_HH__ */
