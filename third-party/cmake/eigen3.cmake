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

add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/eigen3)

set(Eigen3_FOUND TRUE CACHE INTERNAL "" FORCE)
set(EIGEN3_INCLUDE_DIR "${EIGEN3_INCLUDE_DIR};${PYTHON_INCLUDE_DIRS}" CACHE INTERNAL "")
set(EIGEN3_LIBRARIES "${PYTHON_LIBRARIES}" CACHE INTERNAL "")

mask_package_options(EIGEN3)
