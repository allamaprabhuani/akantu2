#===============================================================================
# @file   AkantuMacros.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Wed Feb  9 10:59:42 2011
#
# @brief  Set of macros used by akantu cmake files
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

#===============================================================================
# macro to include optional packages
macro(add_optional_package PACKAGE DESC DEFAULT)
  string(TOUPPER ${PACKAGE} _u_package)
  set(option_name AKANTU_USE_${_u_package})

  option(${option_name} ${DESC} ${DEFAULT})

  if(${option_name})
    if(${PACKAGE} MATCHES BLAS)
      enable_language(Fortran)
    endif()
    find_package(${PACKAGE} REQUIRED)
    if(${_u_package}_FOUND)
      add_definitions(-DAKANTU_USE_${_u_package})
      set(AKANTU_EXTERNAL_LIB_INCLUDE_PATH ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
	${${_u_package}_INCLUDE_PATH}
	)
      set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
	${${_u_package}_LIBRARIES}
	)
      set(AKANTU_${_u_package}_ON ON)
    endif(${_u_package}_FOUND)
  else()
    set(AKANTU_${_u_package}_ON OFF)
  endif(${option_name})
endmacro()

#===============================================================================