if(CMAKE_VERSION VERSION_GREATER 3.8.0)
package_declare(Eigen3 EXTERNAL
  DESCRIPTION "Add Eigen3 dependency to akantu"
  SYSTEM AUTO third-party/cmake/eigen3.cmake
  EXTRA_PACKAGE_OPTIONS ARGS 3.4
  COMPILE_FLAGS CXX -DEIGEN_MAX_CPP_VER=17
  FEATURES_PUBLIC cxx_std_17
  )
else()
package_declare(Eigen3 EXTERNAL
  DESCRIPTION "Add Eigen3 dependency to akantu"
  SYSTEM AUTO third-party/cmake/eigen3.cmake
  EXTRA_PACKAGE_OPTIONS ARGS 3.4
  COMPILE_FLAGS CXX -DEIGEN_MAX_CPP_VER=14
  )
endif()
mark_as_advanced(Eigen3_DIR)
package_add_third_party_script_variable(Eigen3
  EIGEN3_VERSION "3.4.0")
package_add_third_party_script_variable(Eigen3
  EIGEN3_GIT "https://gitlab.com/libeigen/eigen.git")
