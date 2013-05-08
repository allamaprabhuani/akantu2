include(CheckCXXCompilerFlag)
check_cxx_compiler_flag (-std=c++0x HAVE_NEW_STD)

set(AKANTU_CORE_CXX11_FILES
  common/aka_point.hh
  common/aka_bounding_box.hh
  common/aka_bounding_box.cc
  common/aka_geometry.hh
  common/aka_geometry.cc
  model/solid_mechanics/solid_mechanics_model_element.hh
  )


if(HAVE_NEW_STD)
  option(AKANTU_CORE_CXX11 "core CXX11 additions for Akantu" ON)
  add_definitions(-std=c++0x)
else()
  set(AKANTU_CORE_CXX11 OFF CACHE BOOL "core package for Akantu" FORCE)
endif()

mark_as_advanced(AKANTU_CORE_CXX11)
