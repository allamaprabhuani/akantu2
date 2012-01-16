add_optional_package(Mumps "Add Mumps support in akantu" OFF)
set(MUMPS_FILES
  solver/solver_mumps.cc
  solver/solver_mumps.hh
  )


set(MUMPS_DEB_DEPEND
  libmumps-seq-dev
  )
