option(AKANTU_OPTIMIZATION "Use optimization package in Akantu" OFF)

set(AKANTU_OPTIMIZATION_FILES
  common/aka_optimize.hh
  common/aka_optimize.cc
  )

add_external_package_dependencies(optimization nlopt)
mark_as_advanced(AKANTU_OPTIMIZATION)

set(AKANTU_OPTIMIZATION_DOCUMENTATION "
This activates the optimization routines of Akantu. This is currently needed by the 
contact detection algorithms.
" )
