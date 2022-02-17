/**
 * @file   material_von_mises_mazars_non_local.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Dec 17 2020
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
#include "material_damage_non_local.hh"
#include "material_von_mises_mazars.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_VON_MISES_MAZARS_NON_LOCAL_HH_
#define AKANTU_MATERIAL_VON_MISES_MAZARS_NON_LOCAL_HH_

namespace akantu {

/**
 * Material Mazars Non local + Von Mises plasticity
 *
 * parameters in the material files :
 */
template <UInt spatial_dimension>
class MaterialVonMisesMazarsNonLocal
    : public MaterialDamageNonLocal<spatial_dimension,
                                    MaterialVonMisesMazars<spatial_dimension>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MaterialNonLocalParent =
      MaterialDamageNonLocal<spatial_dimension,
                             MaterialVonMisesMazars<spatial_dimension>>;

  MaterialVonMisesMazarsNonLocal(SolidMechanicsModel & model,
                                 const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  void computeNonLocalStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost) override;

  void registerNonLocalVariables() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the ehat per quadrature points to perform the averaging
  InternalField<Real> Ehat;

  InternalField<Real> non_local_variable;
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_VON_MISES_MAZARS_NON_LOCAL_HH_ */
