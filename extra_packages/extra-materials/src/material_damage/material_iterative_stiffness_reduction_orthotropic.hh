/**
 * @file   material_iterative_stiffness_reduction.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
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
#include "material_damage_iterative_orthotropic.hh"
#include "material_iterative_stiffness_reduction.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_ITERATIVE_STIFFNESS_REDUCTION_ORTHOTROPIC_HH__
#define __AKANTU_MATERIAL_ITERATIVE_STIFFNESS_REDUCTION_ORTHOTROPIC_HH__

namespace akantu {

template <UInt spatial_dimension>
class MaterialIterativeStiffnessReductionOrthotropic
    : public MaterialIterativeStiffnessReduction<
          spatial_dimension,
          MaterialDamageIterativeOrthotropic<spatial_dimension>> {
  using parent = MaterialIterativeStiffnessReduction<
      spatial_dimension, MaterialDamageIterativeOrthotropic<spatial_dimension>>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialIterativeStiffnessReductionOrthotropic(SolidMechanicsModel & model,
                                                 const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// update internal field damage and associated principal directions
  UInt updateDamage() override;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_ITERATIVE_STIFFNESS_REDUCTION_ORTHOTROPIC_HH__ */
