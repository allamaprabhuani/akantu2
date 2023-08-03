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

     /Support[Name]
     - /_Nodes_
     - /Connectivities
     -   /_types_
     -   /...
     - /Data
     -   /Data1
     -     /_types_
     -     /...
     /Global Data
   */
  using dumper::HDF5::File;

  if (not h5) {
    auto path = fs::path(directory);
    path /= filename + ".h5";
    h5 = std::make_unique<File>(support, path);
  }

  h5->dump();

  H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
}

} // namespace akantu
