package_declare(HDF5 EXTERNAL
  DESCRIPTION "Add HDF5 support in akantu"
  )

# HDF5_IS_PARALLEL check in parallel if this is activated
package_on_enabled_script(HDF5
  "get_property(_vars DIRECTORY PROPERTY VARIABLES)
  foreach(_var ${_vars})
      message(\"${_var}: ${${_var}}\")
  endforeach()"
)

package_declare_sources(HDF5
  io/new_dumpers/dumper_file_base.hh

  io/new_dumpers/dumper_hdf5.cc
  io/new_dumpers/dumper_hdf5.hh
  io/new_dumpers/dumper_xdmf.cc
  io/new_dumpers/dumper_xdmf.hh

  io/new_dumpers/hdf5_file.hh
  io/new_dumpers/hdf5_entities.hh
  io/new_dumpers/xdmf_file.hh
  io/new_dumpers/xml_helper.hh
  )
