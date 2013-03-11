include(CheckCXXCompilerFlag)
check_cxx_compiler_flag (-std=c++0x HAVE_NEW_STD)

set(AKANTU_CORE_CXX11_FILES
  common/aka_point.hh
  common/aka_bounding_box.hh
  common/aka_bounding_box.cc
  common/aka_geometry.hh
  common/aka_geometry.cc
  )


if(HAVE_NEW_STD)
  set(AKANTU_CORE_CXX11 ON CACHE INTERNAL "core package for Akantu" FORCE)
  add_definitions(-std=c++0x)
endif()
