/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_non_local.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH_
#define AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH_

namespace akantu {

template <Int dim, class MaterialDamageLocal>
class MaterialDamageNonLocal
    : public MaterialNonLocal<dim, MaterialDamageLocal> {
public:
  using MaterialParent = MaterialNonLocal<dim, MaterialDamageLocal>;

  MaterialDamageNonLocal(SolidMechanicsModel & model, const ID & id)
      : MaterialParent(model, id){};

public:
  decltype(auto) getArguments(ElementType el_type, GhostType ghost_type) {
    return MaterialDamageLocal::getArguments(el_type, ghost_type);
  }

protected:
  /* ------------------------------------------------------------------------ */
  virtual void computeNonLocalStress(ElementType type,
                                     GhostType ghost_type = _not_ghost) = 0;

  /* ------------------------------------------------------------------------ */
  void computeNonLocalStresses(GhostType ghost_type) override {
    AKANTU_DEBUG_IN();

    for (auto type : this->getElementFilter().elementTypes(dim, ghost_type)) {
      auto & elem_filter = this->getElementFilter(type, ghost_type);
      if (elem_filter.empty()) {
        continue;
      }

      computeNonLocalStress(type, ghost_type);
    }

    AKANTU_DEBUG_OUT();
  }
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH_ */
