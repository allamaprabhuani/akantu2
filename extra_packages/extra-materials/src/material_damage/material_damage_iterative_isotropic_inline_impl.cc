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
#include "communicator.hh"
#include "material_damage_iterative_isotropic.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
MaterialDamageIterativeIsotropic<spatial_dimension, ElasticParent>::
    MaterialDamageIterativeIsotropic(SolidMechanicsModel & model, const ID & id)
    : parent(model, id) {}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
void MaterialDamageIterativeIsotropic<
    spatial_dimension, ElasticParent>::computeStress(ElementType el_type,
                                                     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  parent::computeStress(el_type, ghost_type);

  auto dam_it = this->damage(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeDamageAndStressOnQuad(sigma, *dam_it);

  ++dam_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->computeNormalizedEquivalentStress(el_type, ghost_type);
  this->norm_max_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
inline void MaterialDamageIterativeIsotropic<spatial_dimension, ElasticParent>::
    computeDamageAndStressOnQuad(Matrix<Real> & sigma, Real & dam) {
  sigma *= 1 - dam;
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
void MaterialDamageIterativeIsotropic<spatial_dimension, ElasticParent>::
    computeTangentModuli(const ElementType & el_type,
                         Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  ElasticParent<spatial_dimension>::computeTangentModuli(
      el_type, tangent_matrix, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

  computeTangentModuliOnQuad(tangent, *dam);
  ++dam;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
void MaterialDamageIterativeIsotropic<spatial_dimension, ElasticParent>::
    computeTangentModuliOnQuad(Matrix<Real> & tangent, Real & dam) {
  tangent *= (1 - dam);
}
/* --------------------------------------------------------------------------
 */

} // namespace akantu
