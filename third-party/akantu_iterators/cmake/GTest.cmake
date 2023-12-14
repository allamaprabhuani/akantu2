if(NOT _GTest_version)
  set(_GTest_version "main")
endif()

include(FetchContent)
FetchContent_Declare(
  googletest
  SYSTEM
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        ${_GTest_version}
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
# set(gtest_build_tests OFF CACHE INTERNAL "" FORCE)
# set(INSTALL_GTEST OFF CACHE INTERNAL "" FORCE)
# set(BUILD_GTEST ON CACHE INTERNAL "" FORCE)
# set(BUILD_GMOCK OFF CACHE INTERNAL "" FORCE)

# set(Python_ADDITIONAL_VERSIONS ${AKANTU_ITERATORS_PYTHON_MAJOR_VERSION})

# add_subdirectory(${_GTest_external_dir}/google-test)

# set_property(TARGET gtest_main PROPERTY CXX_STANDARD 14)
# set_property(TARGET gtest PROPERTY CXX_STANDARD 14)

# add_library(GTest::Main ALIAS gtest_main)
# add_library(GTest::GTest ALIAS gtest)

# mark_as_advanced(INSTALL_GTEST)
# mark_as_advanced_prefix(gtest)
