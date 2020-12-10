/**
 * @file   material_iterative_stiffness_reduction.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Feb 18 15:25:05 2016
 *
 * @brief  Damage material with constant stiffness reduction
 *
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
#include "material_damage_iterative_isotropic.hh"

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_MATERIAL_ITERATIVE_STIFFNESS_REDUCTION_HH_
#define AKANTU_MATERIAL_ITERATIVE_STIFFNESS_REDUCTION_HH_

namespace akantu {

/**
 * Material damage iterative
 *
 * parameters in the material files :
 *   - Gfx
 *   - h
 *   - Sc
 */
/// Proposed by Rots and Invernizzi, 2004: Regularized sequentially linear
// saw-tooth softening model (section 4.2)
template <UInt spatial_dimension, class IsotropicParent =
                                      MaterialDamageIterativeIsotropic<
                                          spatial_dimension, MaterialElastic>>
class MaterialIterativeStiffnessReduction
    : public IsotropicParent {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialIterativeStiffnessReduction(SolidMechanicsModel & model,
                                      const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// init the material
  void initMaterial() override;

  /// update internal field damage
  UInt updateDamage() override;

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the ultimate strain
  InternalField<Real> eps_u;

  /// the initial strength
  InternalField<Real> Sc_init;

  /// the tangent of the tensile stress-strain softening
  InternalField<Real> D;

  /// fracture energy
  Real Gf;

  /// crack_band_width for normalization of fracture energy
  Real crack_band_width;

  /// the reduction constant (denoted by a in the paper of rots)
  Real reduction_constant;
};

template<UInt dim, template<UInt> class ElasticParent>
using MaterialIterativeStiffnessReductionIsotropic = MaterialIterativeStiffnessReduction<dim, MaterialDamageIterativeIsotropic<dim, ElasticParent>>;

} // namespace akantu

#include "material_iterative_stiffness_reduction_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_ITERATIVE_STIFFNESS_REDUCTION_HH__ */
