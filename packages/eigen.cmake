package_declare(Eigen3 EXTERNAL
  DESCRIPTION "Add Eigen3 dependency to akantu"
  SYSTEM AUTO third-party/cmake/eigen3.cmake
  EXTRA_PACKAGE_OPTIONS ARGS 3.4
  )

mark_as_advanced(Eigen3_DIR)
package_add_third_party_script_variable(Eigen3
  EIGEN3_VERSION "3.4.0")
package_add_third_party_script_variable(Eigen3
  EIGEN3_GIT "https://gitlab.com/libeigen/eigen.git")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "12.0.0")
  package_set_compile_flags(Eigen3 CXX "-Wno-array-bounds -Wno-use-after-free -Wno-stringop-overread")
endif()
