/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_mazars.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH_
#define AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH_

namespace akantu {

/**
 * Material Mazars Non local
 *
 * parameters in the material files :
 */
template <Int dim, template <Int> class Parent = MaterialElastic>
class MaterialMazarsNonLocal
    : public MaterialDamageNonLocal<dim, MaterialMazars<dim, Parent>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using parent = MaterialDamageNonLocal<dim, MaterialMazars<dim, Parent>>;

  MaterialMazarsNonLocal(SolidMechanicsModel & model, const ID & id = "");

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
  decltype(auto) getArguments(ElementType el_type, GhostType ghost_type) {
    return zip_replace(parent::getArguments(el_type, ghost_type),
                       "Ehat"_n = make_view(this->Ehat(el_type, ghost_type)));
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the ehat per quadrature points to perform the averaging
  InternalField<Real> Ehat;

  InternalField<Real> non_local_variable;
};

} // namespace akantu

#include "material_mazars_non_local_tmpl.hh"

#endif /* AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH_ */
