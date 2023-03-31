set(_working_dir ${PROJECT_BINARY_DIR}/third-party/src/eigen3-download)
configure_file(${PROJECT_SOURCE_DIR}/third-party/eigen3.cmake.in ${_working_dir}/CMakeLists.txt)

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/third-party/eigen3/CMakeLists.txt)
  message(STATUS "Downloading eigen3")
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/configure-out.log
    ERROR_FILE ${_working_dir}/configure-error.log)

  execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/build-out.log
    ERROR_FILE ${_working_dir}/build-error.log)
endif()

set(CMAKE_BUILD_TYPE Release)
set(BUILD_TESTING OFF CACHE BOOL "Eigen Tests" FORCE)

add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/eigen3)

set(EIGEN3_LIBRARIES Eigen3::Eigen)
package_set_libraries(Eigen3 ${EIGEN3_LIBRARIES})

mask_package_options(EIGEN3)
mask_package_options(EIGEN)
mask_package_options(PASTIX)
mark_as_advanced(
  BUILD_TESTING
  INCLUDE_INSTALL_DIR
  PKGCONFIG_INSTALL_DIR
  PATCH_COMMAND
  QT_QMAKE_EXECUTABLE
  RT_LIBRARY
  CMAKEPACKAGE_INSTALL_DIR
  BLA_VENDOR
  )


include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen3
  REQUIRED_VARS EIGEN3_LIBRARIES
  VERSION_VAR EIGEN3_VERSION)
