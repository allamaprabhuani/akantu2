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
  using dumper::H5File;

  auto path = fs::path(directory);
  path /= filename + ".h5";

  if (not h5) {
    h5 = std::make_unique<H5File>(support, path);
    aka::as_type<H5File>(*h5).create();
  } else {
    aka::as_type<H5File>(*h5).open();
  }

  h5->dump();

  aka::as_type<H5File>(*h5).close();
}

/* -------------------------------------------------------------------------- */
void DumperHDF5::readInternal() {
  using dumper::H5File;

  auto path = fs::path(directory);
  path /= filename + ".h5";

  if (h5) {
    h5.reset();
  }

  h5 = std::make_unique<H5File>(support, path);
  auto & h5file = aka::as_type<H5File>(*h5);

  h5file.open();

  if (not h5file.isAkantuFile()) {
    AKANTU_EXCEPTION(
        "Cannot read this hdf5 file, it was not generated with akantu.");
  }

  h5file.read();

  h5file.close();
}

} // namespace akantu
