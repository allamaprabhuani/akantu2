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

#===============================================================================
# Tests
#===============================================================================

option(BUILD_ALL_TESTS "Build all tests")

macro(REGISTER_TEST test_name source_list)
  add_executable(${test_name} ${source_list})
  target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})
endmacro()

macro(add_mesh MESH_TARGET GEO_FILE DIM ORDER)
  set(_geo_file ${CMAKE_CURRENT_SOURCE_DIR}/${GEO_FILE})
  get_filename_component(_msh_file "${GEO_FILE}" NAME_WE)
  set(_msh_file ${CMAKE_CURRENT_BINARY_DIR}/${_msh_file}.msh)
  if(EXISTS ${_geo_file})
    add_custom_command(
      OUTPUT ${_msh_file}
      DEPENDS ${_geo_file}
      COMMAND gmsh
      ARGS -${DIM} -order ${ORDER} -optimize -o ${_msh_file} ${_geo_file} 2>&1 > /dev/null
      )
    add_custom_target(${MESH_TARGET}
      DEPENDS ${_msh_file})
  else(EXISTS ${_geo_file})
    message(FATAL_ERROR "File ${_geo_file} not found")
  endif(EXISTS ${_geo_file})
endmacro()


macro(ADD_EXAMPLE example_name desc)
  string(TOUPPER ${example_name} upper_name)
  set(AKANTU_UNIT_TESTS "")
  option(BUILD_${upper_name} "${desc}")
  if(BUILD_ALL_TESTS)
    set(BUILD_${upper_name}_OLD ${BUILD_${upper_name}} CACHE INTERNAL "${desc}" FORCE)
    set(BUILD_${upper_name} ON CACHE INTERNAL "${desc}" FORCE)
  else(BUILD_ALL_TESTS)
    if(DEFINED BUILD_${upper_name}_OLD)
      set(BUILD_${upper_name} ${BUILD_${upper_name}_OLD} CACHE BOOL "${desc}" FORCE)
      unset(BUILD_${upper_name}_OLD CACHE)
    endif(DEFINED BUILD_${upper_name}_OLD)
  endif(BUILD_ALL_TESTS)

  if(BUILD_${upper_name})
    add_subdirectory(${example_name})
  endif(BUILD_${upper_name})
endmacro()


add_example(test_facet_extraction "Test mesh utils facet extraction")
add_example(test_fem "Test finite element functionalties")
add_example(test_mesh_io "Test mesh io object")
add_example(test_solid_mechanics_model "Test solid mechanics object")
add_example(test_static_memory "Test static memory")
add_example(test_vector "Test akantu vector")
if(AKANTU_USE_SCOTCH AND SCOTCH_FOUND)
add_example(test_mesh_partitionate "Test mesh partition creation")
endif(AKANTU_USE_SCOTCH AND SCOTCH_FOUND)


#set(AKANTU_UNIT_TESTS
#  test/test_static_memory.cc
#  test/test_mesh_io_msh.cc
#  test/test_vector.cc
#  test/test_fem.cc
#  test/test_solid_mechanics_model.cc
#  test/test_solid_mechanics_model_cube3d.cc
#  test/test_solid_mechanics_model_bar_traction2d.cc
#  test/test_facet_extraction_triangle1.cc
#  test/test_facet_extraction_tetra1.cc
#  test/test_interpolate_line_1.cc
#  test/test_interpolate_triangle_1.cc
#  test/test_interpolate_tetrahedra_1.cc
#  )
#
#if(AKANTU_USE_SCOTCH AND SCOTCH_FOUND)
#  set(AKANTU_UNIT_TESTS
#    ${AKANTU_UNIT_TESTS}
#    test/test_mesh_partitionate.cc
#    )
#endif(AKANTU_USE_SCOTCH AND SCOTCH_FOUND)
#
#foreach(_test ${AKANTU_UNIT_TESTS})
#  get_filename_component(_target "${_test}" NAME_WE)
#  add_executable(${_target} ${_test})
#  target_link_libraries(${_target} akantu ${AKANTU_EXTERNAL_LIBRARIES})
#endforeach(_test)
