/**
 * @file   aka_fwd.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Thu Jan 3 15:15:00 2013
 *
 * @brief  File containing forward declarations in akantu
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef __AKANTU_FWD_HH__
#define __AKANTU_FWD_HH__


__BEGIN_AKANTU__


template<typename T> class Matrix;
template<typename T> class Vector;

template<typename T, bool is_scal = is_scalar<T>::value > class Array;

template <class T> class SpatialGrid;

extern const Array<UInt> empty_filter;

__END_AKANTU__


#endif /* __AKANTU_FWD_HH__ */
