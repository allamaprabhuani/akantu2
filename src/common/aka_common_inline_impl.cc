/**
 * @file   aka_common_inline_impl.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Fri Jun 11 09:48:06 2010
 *
 * @namespace akantu
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
 * @section DESCRIPTION
 *
 * All common things to be included in the projects files
 *
 */

#include <algorithm>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline std::string to_lower(const std::string & str) {
  std::string lstr = str;
  std::transform(lstr.begin(),
		 lstr.end(),
		 lstr.begin(),
		 (int(*)(int))std::tolower);
  return lstr;
}

/* -------------------------------------------------------------------------- */
inline std::string trim(const std::string & to_trim) {
  std::string trimed = to_trim;
  //left trim
  trimed.erase(trimed.begin(),
	       std::find_if(trimed.begin(),
			    trimed.end(),
			    std::not1(std::ptr_fun<int, int>(std::isspace))));
  // right trim
  trimed.erase(std::find_if(trimed.rbegin(),
			    trimed.rend(),
			    std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
	       trimed.end());
  return trimed;
}

__END_AKANTU__
