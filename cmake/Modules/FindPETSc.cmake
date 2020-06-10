# - Try to find PETSc
#  PETSC_FOUND         - system has PETSc
#  PETSC_INCLUDE_DIRS  - the PETSc include directories
#  PETSC_LIBRARIES     - Link these to use PETSc
#  PETSC_VERSION       - Version string (MAJOR.MINOR.SUBMINOR)

if(PETSc_FIND_REQUIRED)
  set(_required REQUIRED)
endif()

find_package(PkgConfig ${_required})

if(NOT PETSc_FIND_REQUIRED AND NOT PKG_CONFIG_FOUND)
  return()
endif()

# reset the search if needed
# get_property(_vars DIRECTORY PROPERTY VARIABLES)
# foreach(_var ${_vars})
#   if (("${_var}" MATCHES "_petsc") OR ("${_var}" MATCHES "PETSC"))
#     message("${_var} -> ${${_var}}")
#     unset(${_var})
#   endif()
# endforeach()

pkg_search_module(_petsc ${_required} PETSc)

if(_petsc_FOUND AND _petsc_VERSION)
  set(PETSC_VERSION ${_petsc_VERSION})
endif()

if(_petsc_FOUND AND (NOT PETSC_LIBRARIES))
  set(_petsc_libs)
  foreach(_lib ${_petsc_LIBRARIES})
    string(TOUPPER "${_lib}" _u_lib)
    find_library(PETSC_LIBRARY_${_u_lib} ${_lib} PATHS ${_petsc_LIBRARY_DIRS})
    list(APPEND _petsc_libs ${PETSC_LIBRARY_${_u_lib}})
    mark_as_advanced(PETSC_LIBRARY_${_u_lib})
  endforeach()

  set(PETSC_LIBRARIES ${_petsc_libs} CACHE FILEPATH "")
  set(PETSC_INCLUDE_DIRS ${_petsc_INCLUDE_DIRS} CACHE PATH "")
  mark_as_advanced(
    PETSC_LIBRARIES
    PETSC_INCLUDE_DIRS
    )

  if(NOT TARGET petsc::petsc)
    add_library(petsc::petsc INTERFACE IMPORTED)
    set_property(TARGET petsc::petsc PROPERTY INTERFACE_LINK_LIBRARIES ${PETSC_LIBRARIES})
    set_property(TARGET petsc::petsc PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS})
  endif()
endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
  REQUIRED_VARS PETSC_LIBRARIES PETSC_INCLUDE_DIRS
  VERSION_VAR PETSC_VERSION)
