add_optional_package(PTScotch "Add PTScotch support in akantu" OFF)
add_optional_package(Scotch "Add Scotch support in akantu" OFF)

set(SCOTCH_FILES
  mesh_utils/mesh_partition/mesh_partition_scotch.cc
  )

set(PTSCOTCH_FILES
  mesh_utils/mesh_partition/mesh_partition_scotch.cc
  )

if(AKANTU_SCOTCH_ON OR AKANTU_PTSCOTCH_ON)
  set(AKANTU_PARTITIONER_ON ON)
else()
  set(AKANTU_PARTITIONER_ON OFF)
endif()

if(AKANTU_PTSCOTCH_ON)
  add_definitions(-DAKANTU_USE_${_u_package})
endif()

