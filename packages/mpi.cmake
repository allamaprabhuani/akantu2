add_optional_package(MPI "Add MPI support in akantu" OFF PREFIX MPI_C MPI DEPENDS SCOTCH)
set(AKANTU_MPI_FILES
  synchronizer/static_communicator_mpi.cc
  synchronizer/static_communicator_mpi_inline_impl.hh
  synchronizer/static_communicator_mpi.hh
  )
