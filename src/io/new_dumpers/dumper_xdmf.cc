/**
 * @file   dumper_xdmf.cc
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
#include "dumper_xdmf.hh"
#include "support.hh"
#include "xdmf_file.hh"
/* -------------------------------------------------------------------------- */
#include <filesystem>
/* -------------------------------------------------------------------------- */

namespace fs = std::filesystem;

namespace akantu {

/* -------------------------------------------------------------------------- */
DumperXdmf::DumperXdmf(dumper::SupportBase & support) : DumperHDF5(support) {}

/* -------------------------------------------------------------------------- */
void DumperXdmf::dumpInternal() {
  if (time_activated) {
    support.addProperty("dump_count", count);
    support.addProperty("time", Real(count));
  }

  DumperHDF5::dumpInternal();
  using dumper::XDMF::File;

  if (not xdmf) {
    auto path = fs::path(directory);
    path /= filename + ".xmf";
    xdmf = std::make_unique<File>(support, path);
  }

  xdmf->dump();
}

} // namespace akantu
