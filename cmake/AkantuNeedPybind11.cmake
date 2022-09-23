if(DEFINED AKANTU_NEED_PYBIND11_LOADED)
  return()
endif()
set(AKANTU_NEED_PYBIND11_LOADED TRUE)


set(PYBIND11_PYTHON_VERSION ${AKANTU_PREFERRED_PYTHON_VERSION} CACHE INTERNAL "")

find_package(pybind11 QUIET)

if (NOT pybind11_FOUND)
  set(PYBIND11_VERSION "v2.4.2")
  set(PYBIND11_GIT "https://github.com/pybind/pybind11.git")

  include(${PROJECT_SOURCE_DIR}/third-party/cmake/pybind11.cmake)
endif()
