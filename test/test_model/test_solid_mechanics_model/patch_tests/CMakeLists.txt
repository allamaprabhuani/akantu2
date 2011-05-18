set(LIST_TYPES
  segment_2
  segment_3
  triangle_3
  triangle_6
  quadrangle_4
  quadrangle_8
  tetrahedron_4
  tetrahedron_10
  hexahedron_8
  )

add_subdirectory(explicit)

if(AKANTU_MUMPS_ON)
  add_subdirectory(implicit)
endif()