#add_optional_package(IOHelper "Add IOHelper support in akantu" ON)

option(AKANTU_USE_IOHELPER "Add IOHelper support in akantu" ON)
mark_as_advanced(AKANTU_USE_IOHELPER)

if(AKANTU_USE_IOHELPER)
  set(IOHELPER_TARGETS_EXPORT AkantuLibraryDepends)
  add_subdirectory(third-party/iohelper)

  list(APPEND AKANTU_EXTERNAL_LIBRARIES iohelper)
  list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third-party/iohelper/src)

  list(APPEND AKANTU_EXPORT_LIST iohelper)

  mark_as_advanced(IOHELPER_TESTS)
  set(AKANTU_IOHELPER ON)
else()
  set(AKANTU_IOHELPER OFF)
endif()


set(AKANTU_IOHELPER_FILES
  io/dumper/dumper_iohelper.hh
  io/dumper/dumper_iohelper.cc
  io/dumper/dumper_iohelper_tmpl.hh
  io/dumper/dumper_paraview.hh
  io/dumper/dumper_paraview.cc
  io/dumper/dumper_iohelper_tmpl_elemental_field.hh
  io/dumper/dumper_iohelper_tmpl_homogenizing_field.hh
  io/dumper/dumper_iohelper_tmpl_material_internal_field.hh
  io/dumper/dumper_iohelper_tmpl_nodal_field.hh
  io/dumper/dumper_iohelper_tmpl_quadrature_points_field.hh
  )
