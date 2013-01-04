option(AKANTU_OPTIMIZATION "Use optimization package in Akantu" OFF)


if (AKANTU_OPTIMIZATION)

  add_optional_external_package(CBLAS "Use CBLAS library" ON)

  message("Change this to an *** INTERNAL *** variable")
  add_optional_external_package(NLopt "Add NLopt optimization support in akantu"  ON)

endif()


set(AKANTU_OPTIMIZATION_FILES
  common/aka_optimize.hh
  common/aka_optimize.cc
  )


