#===============================================================================
# @file   AkantuExtraCompilationProfiles.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Dec 02 2016
# @date last modification: Wed Feb 03 2021
#
# @brief  Compilation profiles
#
#
# @section LICENSE
#
# Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


#Profiling
set(_profiling "-g -ggdb3 -pg -DNDEBUG -DAKANTU_NDEBUG -O2")
set(CMAKE_CXX_FLAGS_PROFILING ${_profiling}
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_C_FLAGS_PROFILING ${_profiling}
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_Fortran_FLAGS_PROFILING ${_profiling}
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "-pg"
  CACHE STRING "Flags used by the linker during profiling builds")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "-pg"
  CACHE STRING "Flags used by the linker during profiling builds")

mark_as_advanced(
  CMAKE_CXX_FLAGS_PROFILING
  CMAKE_C_FLAGS_PROFILING
  CMAKE_Fortran_FLAGS_PROFILING
  CMAKE_EXE_LINKER_FLAGS_PROFILING
  CMAKE_SHARED_LINKER_FLAGS_PROFILING
  )

set(_coverage "-g -ggdb3 -DNDEBUG -DAKANTU_NDEBUG -O2 --coverage")
set(CMAKE_CXX_FLAGS_COVERAGE ${_coverage}
  CACHE STRING "Flags used by the compiler during profiling builds" FORCE)
set(CMAKE_C_FLAGS_COVERAGE ${_coverage}
  CACHE STRING "Flags used by the compiler during profiling builds" FORCE)
set(CMAKE_Fortran_FLAGS_COVERAGE ${_coverage}
  CACHE STRING "Flags used by the compiler during profiling builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE ${_coverage}
  CACHE STRING "Flags used by the compiler during profiling builds" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE ${_coverage}
  CACHE STRING "Flags used by the linker during sanitizing builds" FORCE)

mark_as_advanced(
  CMAKE_CXX_FLAGS_COVERAGE
  CMAKE_C_FLAGS_COVERAGE
  CMAKE_Fortran_FLAGS_COVERAGE
  CMAKE_SHARED_LINKER_FLAGS_COVERAGE
  CMAKE_EXE_LINKER_FLAGS_COVERAGE
  )

# Sanitize the code
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.2") OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(_blacklist " -fsanitize-blacklist=${PROJECT_SOURCE_DIR}/cmake/sanitize-blacklist.txt")
  endif()
  set(_sanitize "-g -ggdb3 -O2 -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer${_blacklist}")

  set(CMAKE_CXX_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_C_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_Fortran_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_EXE_LINKER_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the linker during sanitizing builds")
  set(CMAKE_SHARED_LINKER_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the linker during sanitizing builds")

  mark_as_advanced(
    CMAKE_CXX_FLAGS_SANITIZE
    CMAKE_C_FLAGS_SANITIZE
    CMAKE_Fortran_FLAGS_SANITIZE
    CMAKE_SHARED_LINKER_FLAGS_SANITIZE
    CMAKE_EXE_LINKER_FLAGS_SANITIZE
    )
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(_sanitize "-g -ggdb3 -O2 -fPIE -fsanitize=memory -fsanitize-memory-track-origins -fsanitize-recover=all -fno-omit-frame-pointer -fsanitize-blacklist=${PROJECT_SOURCE_DIR}/cmake/sanitize-blacklist.txt")

  set(CMAKE_CXX_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_C_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_Fortran_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_EXE_LINKER_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the linker during sanitizing builds")
  set(CMAKE_SHARED_LINKER_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the linker during sanitizing builds")

  mark_as_advanced(
    CMAKE_CXX_FLAGS_SANITIZEMEMORY
    CMAKE_C_FLAGS_SANITIZEMEMORY
    CMAKE_Fortran_FLAGS_SANITIZEMEMORY
    CMAKE_SHARED_LINKER_FLAGS_SANITIZEMEMORY
    CMAKE_EXE_LINKER_FLAGS_SANITIZEMEMORY
    )
endif()
