add_optional_package(Mumps "Add Mumps support in akantu" OFF)
set(AKANTU_MUMPS_FILES
  solver/solver_mumps.cc
  solver/solver_mumps.hh
  )


set(AKANTU_MUMPS_DEB_DEPEND
  libmumps-seq-dev
  )
