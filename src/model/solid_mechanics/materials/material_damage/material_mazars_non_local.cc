/**
 * @file   material_mazars_non_local.cc
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Specialization of the material class for the non-local mazars
 * material
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
#include "material_mazars_non_local.hh"

namespace akantu {

template class MaterialMazarsNonLocal<1>;
template class MaterialMazarsNonLocal<2>;
template class MaterialMazarsNonLocal<3>;

template <Int dim> using MaterialMazarsNonLocal_ = MaterialMazarsNonLocal<dim>;

static bool material_is_alocated_mazars_non_local =
    instantiateMaterial<MaterialMazarsNonLocal_>("mazars_non_local");

} // namespace akantu
