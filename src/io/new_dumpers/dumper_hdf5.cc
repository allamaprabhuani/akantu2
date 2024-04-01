/**
 * @file   dumper_hdf5.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Oct 03 2017
 *
 * @brief Dump data in xdmf format
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "dumper_hdf5.hh"
#include "hdf5_file.hh"
#include "support.hh"
#if defined(AKANTU_USE_MPI)
#include "mpi_communicator_data.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <filesystem>
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace fs = std::filesystem;

/* -------------------------------------------------------------------------- */
DumperHDF5::DumperHDF5(dumper::SupportBase & support) : Dumper(support) {}

/* -------------------------------------------------------------------------- */
void DumperHDF5::dumpInternal() {
  /*
     HDF5 structure

     -/steps
     -  /[time]
     -    /[mesh_name]
     -      /topology
     -        /nodes
     -          /position
     -        /connectivities
     -          /[_types_]
     -          /...
     -      /data
     -        /[data_el]
     -          /_types_
     -          /...
     -        /[data_node]
     -        /...
     -      /groups
     -        /[group_name]
     -          /topology
     -        /data
     -          /[data_el]
   */
  using dumper::HDF5::File;

  auto fapl_id = H5Pcreate(H5P_FILE_ACCESS);
#if defined(AKANTU_USE_MPI)
  auto && comm = aka::as_type<MPICommunicatorData>(
                     support.getCommunicator().getCommunicatorData())
                     .getMPICommunicator();
  H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
  /*
   * OPTIONAL: Set collective metadata reads on FAPL to allow
   *           parallel writes to filtered datasets to perform
   *           better at scale. While not strictly necessary,
   *           this is generally recommended.
   */
  H5Pset_all_coll_metadata_ops(fapl_id, true);
#endif

  /*
   * OPTIONAL: Set the latest file format version for HDF5 in
   *           order to gain access to different dataset chunk
   *           index types and better data encoding methods.
   *           While not strictly necessary, this is generally
   *           recommended.
   */
  H5Pset_libver_bounds(fapl_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

  auto path = fs::path(directory);
  path /= filename + ".h5";

  if (not h5) {
    h5 = std::make_unique<File>(support, path, fapl_id);
  } else {
    aka::as_type<File>(h5.get())->open(path, fapl_id);
  }

  h5->dump();

  //  H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

  aka::as_type<File>(h5.get())->close();
  H5Pclose(fapl_id);
}

} // namespace akantu
