#add_optional_package(IOHelper "Add IOHelper support in akantu" ON)

option(AKANTU_USE_IOHELPER "Add IOHelper support in akantu" ON)
mark_as_advanced(AKANTU_USE_IOHELPER)

if(AKANTU_USE_IOHELPER)
  set(IOHELPER_TARGETS_EXPORT AkantuLibraryDepends)
  add_subdirectory(third-party/iohelper)

  list(APPEND AKANTU_EXTERNAL_LIBRARIES iohelper)
  list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${AKANTU_BUILD_IOHELPER_INCLUDE_DIR})

  list(APPEND AKANTU_EXPORT_LIST iohelper)

  mark_as_advanced(IOHELPER_TESTS)
  set(AKANTU_IOHELPER ON)
else()
  set(AKANTU_IOHELPER OFF)
endif()
