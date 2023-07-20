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
#include "material_damage.hh"
#include "material_damage_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef TEST_MATERIAL_DAMAGE_HH_
#define TEST_MATERIAL_DAMAGE_HH_

using namespace akantu;

template <Int dim>
class TestMaterialDamage
    : public MaterialDamageNonLocal<dim, MaterialDamage<dim, MaterialElastic>> {

  using Parent =
      MaterialDamageNonLocal<dim, MaterialDamage<dim, MaterialElastic>>;

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructor */
  /* ------------------------------------------------------------------------ */
public:
  TestMaterialDamage(SolidMechanicsModel & model, const ID & id);
  /* ------------------------------------------------------------------------ */
  /* Methods */
  /* ------------------------------------------------------------------------ */
public:
  void registerNonLocalVariables() final;

  void computeNonLocalStress(ElementType /*element_type*/,
                             GhostType /*ghost_type*/) final{};

  void insertQuadsInNeighborhoods(GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Members */
  /* ------------------------------------------------------------------------ */
private:
  InternalField<Real> & grad_u_nl;
};

#endif /* TEST_MATERIAL_DAMAGE_HH_ */
