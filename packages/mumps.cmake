add_optional_package(Mumps "Add Mumps support in akantu" OFF)
set(MUMPS_FILES
  solver/solver_mumps.cc
  )


set(MUMPS_DEB_DEPEND
  libmumps-seq-dev
  )
