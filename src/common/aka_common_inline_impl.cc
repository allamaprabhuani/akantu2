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

/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING
#include <boost/preprocessor.hpp>

#define AKANTU_BOOST_CASE_MACRO(r,macro,type)	\
  case type : { macro(type); break;}

#define AKANTU_BOOST_ELEMENT_SWITCH(macro);				\
  do {									\
    switch(type) {							\
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO,macro,AKANTU_ELEMENT_TYPE) \
    case _not_defined:							\
    case _max_element_type:  {						\
      AKANTU_DEBUG_ERROR("Wrong type : " << type);			\
      break;								\
    }									\
    }									\
  } while(0)

#define AKANTU_BOOST_LIST_MACRO(r,macro,type)	\
  macro(type)

#define AKANTU_BOOST_ELEMENT_LIST(macro)				\
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO,macro,AKANTU_ELEMENT_TYPE)
