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
#include "material_linear_isotropic_hardening.hh"
#include "material_mazars.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VONMISES_MAZARS_HH__
#define __AKANTU_MATERIAL_VONMISES_MAZARS_HH__

namespace akantu {

template <Int dim>
using MaterialVonMisesMazars =
    MaterialMazars<dim, MaterialLinearIsotropicHardening>;

} // namespace akantu

#endif /* __AKANTU_MATERIAL_VONMISES_MAZARS_HH__ */
