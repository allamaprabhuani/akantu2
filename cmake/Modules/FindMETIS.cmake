#===============================================================================
# Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
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


find_path(METIS_INCLUDE_DIR metis.h
  PATHS "${METIS_DIR}"
  ENV METIS_DIR
  PATH_SUFFIXES include
  )

find_library(METIS_LIBRARY NAMES metis
  PATHS "${METIS_DIR}"
  ENV METIS_DIR
  PATH_SUFFIXES lib
  )

mark_as_advanced(METIS_LIBRARY METIS_INCLUDE_DIR)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
  REQUIRED_VARS
    METIS_LIBRARY
    METIS_INCLUDE_DIR)
