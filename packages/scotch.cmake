add_optional_package(PTScotch "Add PTScotch support in akantu" OFF)
add_optional_package(Scotch "Add Scotch support in akantu" OFF)

if(SCOTCH_INCLUDE_DIR)
  file(STRINGS ${SCOTCH_INCLUDE_DIR}/scotch.h SCOTCH_INCLUDE_CONTENT)
  string(REGEX MATCH "_cplusplus" _match ${SCOTCH_INCLUDE_CONTENT})
  if(_match)
    set(AKANTU_SCOTCH_NO_EXTERN ON)
    list(APPEND AKANTU_DEFINITIONS AKANTU_SCOTCH_NO_EXTERN)
  else()
    set(AKANTU_SCOTCH_NO_EXTERN OFF)
  endif()
endif()


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

