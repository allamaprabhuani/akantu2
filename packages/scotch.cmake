#===============================================================================
# @file   scotch.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for scotch
#
#
# @section LICENSE
#
# Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


package_declare(Scotch EXTERNAL
  DESCRIPTION "Add Scotch support in akantu"
  SYSTEM ON third-party/cmake/scotch.cmake)

package_declare_sources(Scotch
  mesh_utils/mesh_partition/mesh_partition_scotch.cc
  )

package_add_third_party_script_variable(Scotch
  SCOTCH_VERSION "5.1.12b")
package_add_third_party_script_variable(Scotch
  SCOTCH_ARCHIVE_HASH "MD5=e13b49be804755470b159d7052764dc0")
package_add_third_party_script_variable(Scotch
  SCOTCH_ARCHIVE "scotch_${SCOTCH_VERSION}_esmumps.tar.gz")
package_add_third_party_script_variable(Scotch
  SCOTCH_URL "https://gforge.inria.fr/frs/download.php/28978/scotch_${SCOTCH_VERSION}_esmumps.tar.gz")


package_get_option_name(Scotch _opt_name)
package_use_system(Scotch _system)
if(${_opt_name} AND _system)
  include(CheckTypeSize)

  package_get_include_dir(Scotch _include_dir)
  if(_include_dir)
    set(CMAKE_EXTRA_INCLUDE_FILES stdio.h scotch.h)
    set(CMAKE_REQUIRED_INCLUDES ${_include_dir})
    check_type_size("SCOTCH_Num" SCOTCH_NUM)

    if(SCOTCH_NUM AND NOT SCOTCH_NUM EQUAL AKANTU_INTEGER_SIZE)
      math(EXPR _n "${AKANTU_INTEGER_SIZE} * 8")
      message(SEND_ERROR "This version of Scotch cannot be used, it is compiled with the wrong size for SCOTCH_Num."
        "Recompile Scotch with the define -DINTSIZE${_n}. The current scotch integer size is ${SCOTCH_NUM}")
    endif()
  endif()
endif()

package_set_package_system_dependency(Scotch deb libscotch)
package_set_package_system_dependency(Scotch deb-src libscotch-dev)

package_declare_extra_files_to_package(Scotch
  PROJECT
    third-party/cmake/scotch.cmake
    third-party/scotch_5.1.12b.patch
    third-party/scotch_5.1.12b_make.inc.cmake
    cmake/Modules/FindScotch.cmake
  )
