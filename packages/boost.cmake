#===============================================================================
# @file   boost.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for core
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

#set(Boost_USE_STATIC_LIBS        ON)
#set(Boost_USE_MULTITHREADED      ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)

set(AKANTU_BOOST_COMPONENTS
  chrono
#  system
  )


find_package(Boost COMPONENTS ${AKANTU_BOOST_COMPONENTS})

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

foreach(_comp ${AKANTU_BOOST_COMPONENTS})
  string(TOUPPER ${_comp} _u_comp)
  if(Boost_${_u_comp}_FOUND)
    set(AKANTU_${_u_comp}_CHRONO TRUE)
    list(APPEND AKANTU_EXTERNAL_LIBRARIES ${Boost_${_u_comp}_LIBRARY})
  endif()
endforeach()
