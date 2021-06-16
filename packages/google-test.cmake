#===============================================================================
# @file   google-test.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Nov 09 2017
# @date last modification: Fri Mar 16 2018
#
# @brief  package for external dependency to google tests
#
#
# @section LICENSE
#
# Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(gtest EXTERNAL
  DESCRIPTION "Add GTest support for tests"
  SYSTEM AUTO third-party/cmake/gtest.cmake
  EXCLUDE_FROM_ALL
  )

package_get_option_name(gtest _opt_name)
if(AKANTU_TESTS)
  set(${_opt_name} ON CACHE BOOL "Add GTest support for tests (forced)" FORCE)
endif()

package_add_third_party_script_variable(google-test
  GTEST_VERSION "master")
package_add_third_party_script_variable(google-test
  GTEST_GIT "https://github.com/google/googletest.git")
