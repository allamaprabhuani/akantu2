/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <array>
#include <mpi.h>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MPI_TYPE_WRAPPER_HH_
#define AKANTU_MPI_TYPE_WRAPPER_HH_

namespace akantu {

class MPICommunicatorData : public CommunicatorInternalData {
public:
  MPICommunicatorData(const MPI_Comm & comm) {
    if (is_externaly_initialized == -1) {
      MPI_Initialized(&is_externaly_initialized);
    }

    if ((is_externaly_initialized == 0) and (mpi_communicator_instances == 0)) {
      MPI_Init(nullptr, nullptr); // valid according to the spec
    }

    MPI_Comm_create_errhandler(MPICommunicatorData::errorHandler,
                               &error_handler);

    setMPICommunicator(comm);
    ++mpi_communicator_instances;
  }

  MPICommunicatorData(const MPICommunicatorData &) = delete;
  MPICommunicatorData(MPICommunicatorData &&) = delete;
  MPICommunicatorData & operator=(const MPICommunicatorData &) = delete;
  MPICommunicatorData & operator=(MPICommunicatorData &&) = delete;

  ~MPICommunicatorData() override {
    int finalized{0};
    MPI_Finalized(&finalized);

    if ((is_externaly_initialized == 0) and (finalized == 0)) {
      MPI_Comm_set_errhandler(communicator, saved_error_handler);
      MPI_Errhandler_free(&error_handler);

      --mpi_communicator_instances;

      if (mpi_communicator_instances == 0) {
        MPI_Finalize();
      }
    }
  }

  inline void setMPICommunicator(MPI_Comm comm) {
    MPI_Comm_set_errhandler(communicator, saved_error_handler);
    communicator = comm;
    MPI_Comm_get_errhandler(comm, &saved_error_handler);
    MPI_Comm_set_errhandler(comm, error_handler);
  }

  [[nodiscard]] inline int rank() const {
    int prank{};
    MPI_Comm_rank(communicator, &prank);
    return prank;
  }

  [[nodiscard]] inline int size() const {
    int psize{};
    MPI_Comm_size(communicator, &psize);
    return psize;
  }

  [[nodiscard]] inline MPI_Comm getMPICommunicator() const {
    return communicator;
  }

  [[nodiscard]] static int getMaxTag() {
    int flag{};
    int * value{nullptr};
    // not defined on derived intra-communicator
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &value, &flag);
    AKANTU_DEBUG_ASSERT(flag, "No attribute MPI_TAG_UB.");
    return *value;
  }

private:
  MPI_Comm communicator{MPI_COMM_WORLD};
  MPI_Errhandler saved_error_handler{MPI_ERRORS_ARE_FATAL};
  static int
      is_externaly_initialized; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
  static int
      mpi_communicator_instances; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
  /* ------------------------------------------------------------------------ */
  MPI_Errhandler error_handler{};

  static void
  errorHandler(MPI_Comm * /*comm*/,
               int * error_code, // NOLINT(readability-non-const-parameter)
               ...) {
    std::array<char, MPI_MAX_ERROR_STRING> error_string{};
    int str_len{};
    MPI_Error_string(*error_code, error_string.data(), &str_len);
    AKANTU_CUSTOM_EXCEPTION_INFO(debug::CommunicationException(),
                                 "MPI failed with the error code "
                                     << *error_code << ": \""
                                     << error_string.data() << "\"");
  }
};

} // namespace akantu

#endif /* AKANTU_MPI_TYPE_WRAPPER_HH_ */
