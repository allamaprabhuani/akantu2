#===============================================================================
# @file   mumps.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for mumps support
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

package_declare(Mumps EXTERNAL
  DESCRIPTION "Add Mumps support in akantu"
  SYSTEM ON third-party/cmake/mumps.cmake
  )

package_declare_sources(Mumps
  solver/solver_mumps.cc
  solver/solver_mumps.hh
  )

package_get_option_name(parallel _par_option)
if(${_par_option})
  package_set_find_package_extra_options(Mumps ARGS COMPONENTS "parallel")
  package_add_third_party_script_variable(Mumps MUMPS_TYPE "par")

  package_set_package_system_dependency(Mumps deb libmumps)
  package_set_package_system_dependency(Mumps deb-src libmumps-dev)
else()
  package_set_find_package_extra_options(Mumps ARGS COMPONENTS "sequential")
  package_add_third_party_script_variable(Mumps MUMPS_TYPE "seq")

  package_set_package_system_dependency(Mumps deb libmumps-seq)
  package_set_package_system_dependency(Mumps deb-src libmumps-seq-dev)
endif()

package_use_system(Mumps _use_system)
if(NOT _use_system)
  enable_language(Fortran)

  set(AKANTU_USE_MUMPS_VERSION "4.10.0" CACHE STRING "Default Mumps version to compile")
  mark_as_advanced(AKANTU_USE_MUMPS_VERSION)
  set_property(CACHE AKANTU_USE_MUMPS_VERSION PROPERTY STRINGS "4.9.2" "4.10.0" "5.0.0")

  package_get_option_name(MPI _mpi_option)
  if(${_mpi_option})
    package_add_dependencies(Mumps ScaLAPACK MPI)
  endif()

  package_add_dependencies(Mumps Scotch BLAS)
endif()

package_declare_documentation(Mumps
  "This package enables the \\href{http://mumps.enseeiht.fr/}{MUMPS} parallel direct solver for sparce matrices."
  "This is necessary to solve static or implicit problems."
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  ""
  "\\begin{command}"
  "  > sudo apt-get install libmumps-seq-dev # for sequential"
  "  > sudo apt-get install libmumps-dev     # for parallel"
  "\\end{command}"
  ""
  "Under Mac OS X the installation requires the following steps:"
  "\\begin{command}"
  "  > sudo port install mumps"
  "\\end{command}"
  ""
  "If you activate the advanced option AKANTU\\_USE\\_THIRD\\_PARTY\\_MUMPS the make system of akantu can automatically compile MUMPS. For this you will have to download MUMPS from \\url{http://mumps.enseeiht.fr/} or \\url{http://graal.ens-lyon.fr/MUMPS} and place it in \\shellcode{<akantu source>/third-party}"
  )

package_declare_extra_files_to_package(MUMPS
  PROJECT
    third-party/MUMPS_4.10.0_make.inc.cmake
    third-party/MUMPS_5.0.0.patch
    third-party/MUMPS_4.10.0.patch
    third-party/MUMPS_4.9.2_make.inc.cmake
    third-party/cmake/mumps.cmake
    cmake/Modules/FindMumps.cmake
  )
