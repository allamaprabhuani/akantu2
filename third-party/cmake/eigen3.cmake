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

set_property(TARGET eigen APPEND
  PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third-party/eigen3>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

set(Eigen3_FOUND TRUE CACHE INTERNAL "" FORCE)
set(EIGEN3_LIBRARIES eigen CACHE INTERNAL "")

mark_as_advanced(BUILD_TESTING)
mask_package_options(EIGEN3)
mask_package_options(EIGEN)
