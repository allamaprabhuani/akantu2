/**
 * @file   aka_common_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Dec 01 12:54:29 2011
 *
 * @brief  inline implementations of common akantu type descriptions
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
//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementType type)
{
  switch(type)
    {
    case _segment_2        : stream << "_segment_2"       ; break;
    case _segment_3        : stream << "_segment_3"       ; break;
    case _triangle_3       : stream << "_triangle_3"      ; break;
    case _triangle_6       : stream << "_triangle_6"      ; break;
    case _tetrahedron_4    : stream << "_tetrahedron_4"   ; break;
    case _tetrahedron_10   : stream << "_tetrahedron_10"  ; break;
    case _quadrangle_4     : stream << "_quadrangle_4"    ; break;
    case _quadrangle_8     : stream << "_quadrangle_8"    ; break;
    case _hexahedron_8     : stream << "_hexahedron_8"    ; break;
    case _bernoulli_beam_2 : stream << "_bernoulli_beam_2"; break;
#if defined(AKANTU_COHESIVE_ELEMENT)
    case _cohesive_2d_4    : stream << "_cohesive_2d_4"   ; break;
    case _cohesive_2d_6    : stream << "_cohesive_2d_6"   ; break;
#endif
    case _not_defined      : stream << "_not_defined"     ; break;
    case _max_element_type : stream << "ElementType(" << (int) type << ")"; break;
    case _point            : stream << "point"; break;
    }
  return stream;
}

/// standard output stream operator for GhostType
inline std::ostream & operator <<(std::ostream & stream, GhostType type)
{
  switch(type)
    {
    case _not_ghost : stream << "not_ghost"; break;
    case _ghost     : stream << "ghost"    ; break;
    case _casper    : stream << "Casper the friendly ghost"; break;
    }
  return stream;
}

/// standard output stream operator for SynchronizationTag
inline std::ostream & operator <<(std::ostream & stream, SynchronizationTag type)
{
  switch(type)
    {
    case _gst_smm_mass                 : stream << "_gst_smm_mass"                ; break;
    case _gst_smm_for_strain	       : stream << "_gst_smm_for_strain"	  ; break;
    case _gst_smm_boundary	       : stream << "_gst_smm_boundary"	      	  ; break;
    case _gst_smm_uv		       : stream << "_gst_smm_uv"		  ; break;
    case _gst_smm_res		       : stream << "_gst_smm_res"		  ; break;
    case _gst_smm_init_mat	       : stream << "_gst_smm_init_mat"	      	  ; break;
    case _gst_smm_stress	       : stream << "_gst_smm_stress"	      	  ; break;
    case _gst_htm_capacity	       : stream << "_gst_htm_capacity" 	      	  ; break;
    case _gst_htm_temperature	       : stream << "_gst_htm_temperature" 	  ; break;
    case _gst_htm_gradient_temperature : stream << "_gst_htm_gradient_temperature"; break;
    case _gst_mnl_for_average	       : stream << "_gst_mnl_for_average"	  ; break;
    case _gst_mnl_weight               : stream << "_gst_mnl_weight"       	  ; break;
    case _gst_test                     : stream << "_gst_test"                    ; break;
    }
  return stream;
}

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
