option(AKANTU_OPTIMIZATION "Use optimization package in Akantu" OFF)

set(AKANTU_OPTIMIZATION_FILES
  common/aka_optimize.hh
  common/aka_optimize.cc
  )

add_optional_external_package(CBLAS "Use CBLAS library" OFF)
add_optional_external_package(NLopt "Add NLopt optimization support in akantu"  OFF)

add_external_package_dependencies(optimization cblas)

mark_as_advanced(AKANTU_OPTIMIZATION)
mark_as_advanced(AKANTU_USE_CBLAS)
mark_as_advanced(AKANTU_USE_NLOPT)

if(AKANTU_OPTIMIZATION)
  list(APPEND AKANTU_BOOST_COMPONENTS
    chrono
    system
    )
endif()
