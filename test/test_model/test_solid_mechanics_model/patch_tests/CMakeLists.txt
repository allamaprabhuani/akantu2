add_subdirectory(explicit)

if(AKANTU_MUMPS_ON)
  add_subdirectory(implicit)
endif()