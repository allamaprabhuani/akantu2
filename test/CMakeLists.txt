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
# List of tests
#===============================================================================
add_akantu_test(test_surface_extraction "Test mesh utils surface extraction")
add_akantu_test(test_facet_extraction "Test mesh utils facet extraction")
add_akantu_test(test_fem "Test finite element functionalties")
add_akantu_test(test_mesh_io "Test mesh io object")
add_akantu_test(test_solid_mechanics_model "Test solid mechanics object")
add_akantu_test(test_static_memory "Test static memory")
add_akantu_test(test_vector "Test akantu vector")

if(AKANTU_SCOTCH_ON)
  add_akantu_test(test_mesh_partitionate "Test mesh partition creation")
endif(AKANTU_SCOTCH_ON)

if(AKANTU_MPI_ON)
  add_akantu_test(test_synchronizer "Test synchronizers")
endif(AKANTU_MPI_ON)
