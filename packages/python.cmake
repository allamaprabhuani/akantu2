#===============================================================================
# @file   python_interpreter.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Oct 31 2018
# @date last modification: Wed Oct 31 2018
#
# @brief  packages description for the external dependency to the python interpreted
#
#
# @section LICENSE
#
# Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
if(FALSE AND CMAKE_VERSION VERSION_GREATER 3.21)
  if (PYTHON_EXECUTABLE)
    get_filename_component(_python_bin_dir "${PYTHON_EXECUTABLE}" DIRECTORY)
    get_filename_component(_python_root_dir "${_python_bin_dir}" DIRECTORY)
    set(Python_ROOT_DIR "${_python_bin_dir}" CACHE INTERNAL "")
    message(STATUS "Python hint: ${Python_ROOT_DIR}")
  endif()


  package_declare(Python EXTERNAL DESCRIPTION "Akantu's python dependency"
    EXTRA_PACKAGE_OPTIONS ARGS
    COMPONENTS Interpreter;Development
    )

  package_on_enabled_script(Python
    "set(Python_SITELIB \${Python_SITELIB} CACHE INTERNAL \"\")
set(Python_EXECUTABLE \${Python_EXECUTABLE} CACHE INTERNAL \"\")
set(Python_INCLUDE_DIRS \${Python_INCLUDE_DIRS} CACHE INTERNAL \"\")
set(Python_VERSION_MAJOR \${Python_VERSION_MAJOR} CACHE INTERNAL \"\")
set(Python_VERSION_MINOR \${Python_VERSION_MINOR} CACHE INTERNAL \"\")
")
else()
  if (PYTHON_LIBRARY MATCHES "\.a$")
    set(PYTHON_LIBRARY NOTFOUND CACHE INTERNAL "")
  endif()
  package_declare(Python ADVANCED META DESCRIPTION "Akantu's python dependency"
    DEPENDS PythonLibsNew PythonInterp
    )

  package_declare(PythonInterp EXTERNAL
    DESCRIPTION "Akantu's python Interpreter"
    )

  package_on_enabled_script(PythonInterp
    "set(Python_EXECUTABLE \${PYTHON_EXECUTABLE} CACHE INTERNAL \"\")")

  package_declare(PythonLibsNew EXTERNAL
    DESCRIPTION "Akantu's python Development"
    )

  package_on_enabled_script(PythonLibsNew
    "set(Python_INCLUDE_DIRS \${PYTHON_INCLUDE_DIRS} CACHE INTERNAL \"\")
set(Python_VERSION_MAJOR \${PYTHON_VERSION_MAJOR} CACHE INTERNAL \"\")
set(Python_VERSION_MINOR \${PYTHON_VERSION_MINOR} CACHE INTERNAL \"\")
")

  # package_declare(Numpy EXTERNAL
  #   DESCRIPTION "Akantu's python Numpy"
  #   )
endif()
