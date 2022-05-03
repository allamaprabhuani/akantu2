/**
 * @file   material_damage_explicit.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 *
 * @brief  Specialization of the class material damage explicit to damage all
 * the quadrature point above damaging criteria at once. It may inherit both
 * from material elastic and viscoelastic. Max principal stress criterion is
 * used as a failure criterion. Can be used ONLY IN EXPLICIT DYNAMIC SIMULATION.
 * Implementation is based on the one from Code_Aster.
 * https://www.code-aster.org/V2/doc/v11/fr/man_r/r5/r5.03.18.pdf
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_DAMAGE_EXPLICIT_HH__
#define __AKANTU_MATERIAL_DAMAGE_EXPLICIT_HH__

namespace akantu {

/**
 * Material damage explicit
 *
 * parameters in the material files :
 *   - Gf
 *   - h
 *   - Sc
 */
template <UInt dim, template <UInt> class ElasticParent = MaterialElastic>
class MaterialDamageExplicit : public MaterialDamage<dim, ElasticParent> {
  using parent = MaterialDamage<dim, ElasticParent>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDamageExplicit(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// init the material
  void initMaterial() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

protected:
  inline void computeDamageAndStressOnQuad(Matrix<Real> & grad_u,
                                           Matrix<Real> & sigma, Real & dam,
                                           const Real & Et, const Real & Sc);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  using voigt_h = VoigtHelper<dim>;

  /// resistance to damage
  RandomInternalField<Real> Sc;

  /// fracture energy
  Real Gf;

  /// crack_band_width for normalization of fracture energy
  Real crack_band_width;

  /// maximum damage. introduced to avoid zero stiffness
  Real max_damage;

  /// the softening slope
  InternalField<Real> Et;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_damage_explicit_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_DAMAGE_EXPLICIT_HH__ */
