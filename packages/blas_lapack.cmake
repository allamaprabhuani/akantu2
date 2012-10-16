add_optional_package(BLAS "Use BLAS for arithmetic operations" OFF LANGUAGE Fortran)
if(BLAS_mkl_core_LIBRARY)
  set(AKANTU_USE_BLAS_MKL)
endif()

add_optional_package(LAPACK "Use LAPACK for arithmetic operations" OFF LANGUAGE Fortran)
