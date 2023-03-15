/**
 * @file   material_anisotropic_damage.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Jul 03 2019
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Base class for anisotropic damage materials
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_anisotropic_damage.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */

static bool material_is_alocated_anisotropic_damage [[gnu::unused]] =
    MaterialFactory::getInstance().registerAllocator(
        "anisotropic_damage",
        [](Int dim, const ID & option, SolidMechanicsModel & model,
           const ID & id) -> std::unique_ptr<Material> {
          return tuple_dispatch<AllSpatialDimensions>(
              [&](auto && _) -> std::unique_ptr<Material> {
                constexpr auto dim_ = aka::decay_v<decltype(_)>;

                if (option.empty() or option == "mazars") {
                  return std::make_unique<MaterialAnisotropicDamage<
                      dim_, EquivalentStrainMazars, DamageThresholdTan>>(model,
                                                                         id);
                }
                if (option == "mazars-drucker-prager") {
                  return std::make_unique<MaterialAnisotropicDamage<
                      dim_, EquivalentStrainMazarsDruckerPrager,
                      DamageThresholdTan>>(model, id);
                }
                AKANTU_EXCEPTION("The option "
                                 << option << " is not valid for the material "
                                 << id);
              },
              dim);
        });

} // namespace akantu
