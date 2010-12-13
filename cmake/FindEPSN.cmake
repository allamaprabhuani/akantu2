#===============================================================================
# @file   FindEPSN.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Tue Aug  25 16:53:57 2010
#
# @brief  The find_package file for EPSN
#
# @section LICENSE
#
# <insert license here>
#
#===============================================================================

#===============================================================================
find_path(EPSN_DIR EPSNConfig.cmake
  PATHS $ENV{EPSN_TOP}
  )


if(EPSN_DIR)
  include(${EPSN_DIR}/EPSNConfig.cmake)
  set(EPSN_LIB_PATH ${EPSN_DIR}/lib
    ${EPSN_LIBRARIES_DIR}
    )
  find_library(EPSN_COMMON_LIBRARY epsn_common
    PATHS ${EPSN_LIB_PATH}
    )
  find_library(EPSN_SIMULATION_LIBRARY epsn_simulation
    PATHS ${EPSN_LIB_PATH}
    )
  include(${EPSN_DIR}/EPSNLibraryDepends.cmake)

  set(EPSN_LIBRARIES ${EPSN_SIMULATION_LIBRARY} ${EPSN_COMMON_LIBRARY})
endif(EPSN_DIR)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EPSN DEFAULT_MSG
  EPSN_LIBRARIES EPSN_INCLUDE_DIR)
