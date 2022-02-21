#===============================================================================
# @file   blas.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Oct 16 2012
# @date last modification: Thu May 11 2017
#
# @brief  package description for blas support
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


package_declare(BLAS EXTERNAL
  DESCRIPTION "Use BLAS for arithmetic operations"
  EXTRA_PACKAGE_OPTIONS LANGUAGE Fortran)

set(_default_blas $ENV{BLA_VENDOR})
if(NOT _default_blas)
  set(_default_blas All)
endif()

set(AKANTU_USE_BLAS_VENDOR "${_default_blas}" CACHE STRING "Version of blas to use")
mark_as_advanced(AKANTU_USE_BLAS_VENDOR)
set_property(CACHE AKANTU_USE_BLAS_VENDOR PROPERTY STRINGS
  All
  ACML
  ACML_GPU
  ACML_MP
  ATLAS
  Apple
  CXML
  DXML
  Generic
  Goto
  IBMESSL
  Intel
  Intel10_32
  Intel10_64lp
  Intel10_64lp_seq
  NAS
  OpenBLAS
  PhiPACK
  SCSL
  SGIMATH
  SunPerf
  )

set(ENV{BLA_VENDOR} ${AKANTU_USE_BLAS_VENDOR})

if(BLAS_mkl_core_LIBRARY)
  set(AKANTU_USE_BLAS_MKL CACHE INTERNAL "" FORCE)
endif()

package_set_package_system_dependency(BLAS deb libblas3)
package_set_package_system_dependency(BLAS deb-src libblas3)
