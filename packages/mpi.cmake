add_optional_package(MPI "Add MPI support in akantu" OFF)

set(MPI_FILES
  synchronizer/static_communicator_mpi.cc
  synchronizer/static_communicator_mpi_inline_impl.cc
  synchronizer/static_communicator_mpi.hh
  )

set(MPI_DEPENDS
  SCOTCH
  )