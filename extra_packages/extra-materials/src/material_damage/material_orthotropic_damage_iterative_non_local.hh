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
#include "material_damage_non_local.hh"
#include "material_orthotropic_damage_iterative.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_NON_LOCAL_HH_
#define AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_NON_LOCAL_HH_

namespace akantu {

/**
 * Material Damage Iterative Non local
 *
 * parameters in the material files :
 */
template <Int spatial_dimension>
class MaterialOrthotropicDamageIterativeNonLocal
    : public MaterialDamageNonLocal<
          spatial_dimension,
          MaterialOrthotropicDamageIterative<spatial_dimension>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialDamageNonLocal<
      spatial_dimension, MaterialOrthotropicDamageIterative<spatial_dimension>>
      MaterialOrthotropicDamageIterativeNonLocalParent;
  MaterialOrthotropicDamageIterativeNonLocal(SolidMechanicsModel & model,
                                             const ID & id = "");

  virtual ~MaterialOrthotropicDamageIterativeNonLocal(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

protected:
  void computeStress(ElementType type, GhostType ghost_type);

  void computeNonLocalStress(ElementType type,
                             GhostType ghost_type = _not_ghost);

  /// associate the non-local variables of the material to their neighborhoods
  virtual void nonLocalVariableToNeighborhood();

private:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  InternalField<Real> grad_u_nl;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_orthotropic_damage_iterative_non_local_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_NON_LOCAL_HH_ */
