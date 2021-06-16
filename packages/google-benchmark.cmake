#===============================================================================
# @file   google-benchmark.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Jun 19 2018
# @date last modification: Tue Oct 09 2018
#
# @brief  package for external dependency to google benchmarks
#
#
# @section LICENSE
#
# Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(gbenchmark EXTERNAL
  DESCRIPTION "Add Google Benchmark support"
  SYSTEM AUTO third-party/cmake/gbenchmark.cmake
  EXCLUDE_FROM_ALL
  )

package_get_option_name(gbenchmark _opt_name)
package_add_third_party_script_variable(gbenchmark
  GBENCHMARK_VERSION "master")
package_add_third_party_script_variable(gbenchmark
  GBENCHMARK_GIT "https://github.com/google/benchmark.git")
