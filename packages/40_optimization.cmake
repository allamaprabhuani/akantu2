option(AKANTU_OPTIMIZATION "Use optimization package in Akantu" OFF)

set(AKANTU_OPTIMIZATION_FILES
  common/aka_optimize.hh
  common/aka_optimize.cc
  )

add_internal_package_dependencies(optimization nlopt)
mark_as_advanced(AKANTU_OPTIMIZATION)

