#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

INCLUDE(${AKANTU_CMAKE_DIR}/AkantuTestAndExamples.cmake)

#===============================================================================
# List of examples
#===============================================================================
add_example(2d_compression "Compression example")

if(AKANTU_SCOTCH_ON AND AKANTU_MPI_ON)
  add_example(3d_cube_compression "Cube compression example")
endif(AKANTU_SCOTCH_ON AND AKANTU_MPI_ON)

if(AKANTU_EPSN_ON)
  add_example(epsn "Example of steering akantu with EPSN")
endif(AKANTU_EPSN_ON)