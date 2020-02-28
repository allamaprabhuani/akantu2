package_declare(Eigen3 EXTERNAL
  DESCRIPTION "Add Eigen3 dependency to akantu"
  EXTRA_PACKAGE_OPTIONS ARGS 3.3
  COMPILE_FLAGS CXX -DEIGEN_MAX_CPP_VER=14
  )

mark_as_advanced(Eigen3_DIR)
