#add_optional_package(IOHelper "Add IOHelper support in akantu" ON)

option(AKANTU_USE_IOHELPER "Add IOHelper support in akantu" ON)
mark_as_advanced(AKANTU_USE_IOHELPER)

if(AKANTU_USE_IOHELPER)
  set(IOHELPER_TARGETS_EXPORT AkantuLibraryDepends)
  add_subdirectory(third-party/iohelper)
  list(APPEND AKANTU_EXTERNAL_LIBRARIES iohelper)
  mark_as_advanced(IOHELPER_TESTS)
endif()


