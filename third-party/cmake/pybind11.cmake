set(_working_dir ${PROJECT_BINARY_DIR}/third-party/src/pybind11-download)
configure_file(${PROJECT_SOURCE_DIR}/third-party/pybind11.cmake.in ${_working_dir}/CMakeLists.txt)

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/third-party/pybind11/CMakeLists.txt)
  message(STATUS "Downloading pybind11")
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/configure-out.log
    ERROR_FILE ${_working_dir}/configure-error.log)

  execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY ${_working_dir}
    OUTPUT_FILE ${_working_dir}/build-out.log
    ERROR_FILE ${_working_dir}/build-error.log)
endif()

set(PYBIND11_PYTHON_VERSION ${AKANTU_PREFERRED_PYTHON_VERSION})
add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/pybind11)

set(pybind11_FOUND TRUE CACHE INTERNAL "" FORCE)
set(PYBIND11_INCLUDE_DIR "${PYBIND11_INCLUDE_DIR};${PYTHON_INCLUDE_DIRS}" CACHE INTERNAL "")
set(PYBIND11_LIBRARIES "${PYTHON_LIBRARIES}" CACHE INTERNAL "")

mask_package_options(PYBIND11)
mark_as_advanced(
  USE_PYTHON_INCLUDE_DIR
  )

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(pybind11
  REQUIRED_VARS PYBIND11_LIBRARIES PYBIND11_INCLUDE_DIR
  VERSION_VAR PYBIND11_VERSION)
