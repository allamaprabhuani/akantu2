/**
 * @file   material_damage_iterative.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date   Thu Feb 18 15:25:05 2016
 *
 * @brief  Damage material with constant stiffness reduction
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
#include "material_damage_iterative.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_DAMAGE_ITERATIVE_ISOTROPIC_HH__
#define __AKANTU_MATERIAL_DAMAGE_ITERATIVE_ISOTROPIC_HH__

namespace akantu {

/**
 * Material damage iterative isotropic
 *
 */

template <UInt spatial_dimension,
          template <UInt> class ElasticParent = MaterialElastic>
class MaterialDamageIterativeIsotropic
    : public MaterialDamageIterative<spatial_dimension, ElasticParent> {
  using parent = MaterialDamageIterative<spatial_dimension, ElasticParent>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDamageIterativeIsotropic(SolidMechanicsModel & model,
                                   const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeStress(ElementType el_type, GhostType ghost_type) override;

  void computeTangentModuli(ElementType el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type) override;

protected:
  // simple multiplication by the (1 - damage)
  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma, Real & dam);
  // simple multiplication by the (1 - damage)
  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent, Real & dam);
};

} // namespace akantu

#include "material_damage_iterative_isotropic_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE_ISOTROPIC_HH__ */
