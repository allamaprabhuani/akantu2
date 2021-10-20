/**
 * @file   aka_fortran_mangling.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Feb 08 2021
 * @date last modification:  Mon Feb 08 2021
 *
 * @brief  Copy of the auto generated aka_fortran_mangling.hh
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

#ifndef AKA_FC_HEADER_INCLUDED
#define AKA_FC_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define AKA_FC_GLOBAL(name, NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define AKA_FC_GLOBAL_(name, NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define AKA_FC_MODULE(mod_name, name, mod_NAME, NAME) __##mod_name##_MOD_##name

/* Mangling for Fortran module symbols with underscores. */
#define AKA_FC_MODULE_(mod_name, name, mod_NAME, NAME) __##mod_name##_MOD_##name

#endif
