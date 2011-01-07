#===============================================================================
# @file   FindScotch.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Tue Aug  25 16:53:57 2010
#
# @brief  The find_package file for Scotch
#
# @section LICENSE
#
# <insert license here>
#
#===============================================================================

#===============================================================================
#if(SCOTCH_DIR)
#  set(SCOTCH_LIBRARY "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
#endif(SCOTCH_DIR)

find_library(SCOTCH_LIBRARY scotch
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(SCOTCH_LIBRARY_ERR scotcherr
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_path(SCOTCH_INCLUDE_PATH scotch.h
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES include scotch src/libscotch include/scotch
  )

#===============================================================================
mark_as_advanced(SCOTCH_LIBRARY)
mark_as_advanced(SCOTCH_LIBRARY_ERR)
mark_as_advanced(SCOTCH_INCLUDE_PATH)

set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY} ${SCOTCH_LIBRARY_ERR})
set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES_ALL} CACHE INTERNAL "Libraries for scotch" FORCE)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCOTCH DEFAULT_MSG
  SCOTCH_LIBRARY SCOTCH_LIBRARY_ERR SCOTCH_INCLUDE_PATH)

#===============================================================================
if(NOT SCOTCH_FOUND)
  set(SCOTCH_DIR "" CACHE PATH "Location of Scotch library.")
endif(NOT SCOTCH_FOUND)
