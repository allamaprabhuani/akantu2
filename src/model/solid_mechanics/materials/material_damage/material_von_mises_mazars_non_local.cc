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
#include "material_von_mises_mazars_non_local.hh"

namespace akantu {

template class MaterialMazarsNonLocal<1, MaterialLinearIsotropicHardening>;
template class MaterialMazarsNonLocal<2, MaterialLinearIsotropicHardening>;
template class MaterialMazarsNonLocal<3, MaterialLinearIsotropicHardening>;

const bool material_is_allocated_plastic_mazars_non_local [[maybe_unused]] =
    instantiateMaterial<MaterialVonMisesMazarsNonLocal>(
        "plastic_mazars_non_local");

} // namespace akantu
