include(FetchContent)
FetchContent_Declare(
  googletest
  SYSTEM
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        ${_GTest_version}
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


# set(_working_dir ${PROJECT_BINARY_DIR}/third-party/src/gtest-download)

# configure_file(${PROJECT_SOURCE_DIR}/third-party/gtest.cmake.in ${_working_dir}/CMakeLists.txt)

# set(GTEST_ROOT ${PROJECT_BINARY_DIR}/third-party)

# if(NOT EXISTS ${PROJECT_SOURCE_DIR}/third-party/google-test/CMakeLists.txt)
#   message(STATUS "Downloading googletest")
#   execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#     RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
#     OUTPUT_FILE ${_working_dir}/configure-out.log
#     ERROR_FILE ${_working_dir}/configure-error.log)

#   if(result)
#     message(SEND_ERROR "CMake step for googletest failed: ${result}")
#     return()
#   endif()

#   execute_process(COMMAND ${CMAKE_COMMAND} --build .
#     RESULT_VARIABLE result WORKING_DIRECTORY ${_working_dir}
#     OUTPUT_FILE ${_working_dir}/build-out.log
#     ERROR_FILE ${_working_dir}/build-error.log)

#   if(result)
#     message(SEND_ERROR "Downloading googletest failed: ${result}")
#     return()
#   endif()
# endif()

# set(gtest_build_tests OFF CACHE INTERNAL "" FORCE)
# set(gtest_build_samples OFF CACHE INTERNAL "" FORCE)
# set(BUILD_GTEST ON CACHE INTERNAL "" FORCE)
# set(BUILD_GMOCK OFF CACHE INTERNAL "" FORCE)

# set(Python_ADDITIONAL_VERSIONS ${AKANTU_PREFERRED_PYTHON_VERSION})
# add_subdirectory(third-party/google-test)

# set_property(TARGET gtest_main PROPERTY CXX_STANDARD ${AKANTU_CXX_STANDARD})
# set_property(TARGET gtest PROPERTY CXX_STANDARD ${AKANTU_CXX_STANDARD})

# if (NOT TARGET GTest::gtest_main)
#   add_library(GTest::gtest_main ALIAS gtest_main)
# endif()

# set(gtest_FOUND TRUE CACHE INTERNAL "" FORCE)
set(GTEST_LIBRARIES GTest::gtest_main CACHE INTERNAL "" FORCE)

# mark_as_advanced(
#   INSTALL_GTEST
#   GTEST_HAS_ABSL
#   )
# mask_package_options(gtest)
