/**
 * @file   dumper_iohelper.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 26 2012
 * @date last modification: Thu Sep 17 2015
 *
 * @brief  implementation of Dumper
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "dumper.hh"
#include "dumper_field.hh"
#include "dumper_variable.hh"
/* -------------------------------------------------------------------------- */
#include <filesystem>
/* -------------------------------------------------------------------------- */

namespace akantu {

using dumper::SupportBase;

/* -------------------------------------------------------------------------- */
Dumper::Dumper(dumper::SupportBase & support) : support(support) {
  if (support.isDistributed()) {
    auto && communicator = support.getCommunicator();
    prank = communicator.whoAmI();
    psize = communicator.getNbProc();
  }
}

/* -------------------------------------------------------------------------- */
Dumper::~Dumper() = default;

/* -------------------------------------------------------------------------- */
void Dumper::dump() {
  namespace fs = std::filesystem;
  auto path = fs::path(directory);
  if (not fs::exists(path)) {
    fs::create_directories(path);
  }

  for (auto && [_, field] : support.getFields()) {
    field->update();
  }

  for (auto && [_, sub_support] : support.getSubSupports()) {
    for (auto && [_, field] : sub_support->getFields()) {
      field->update();
    }
  }

  dumpInternal();
  if (time_activated) {
    ++count;
    time += time_step;
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::read() {
  namespace fs = std::filesystem;

  // OAreadInternal();

  for (auto && [_, field] : support.getFields()) {
    field->back_propagate();
  }

  for (auto && [_, sub_support] : support.getSubSupports()) {
    for (auto && [_, field] : sub_support->getFields()) {
      field->back_propagate();
    }
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::setDirectory(const std::string & directory) {
  this->directory = directory;
}

/* -------------------------------------------------------------------------- */
void Dumper::setBaseName(const std::string & basename) {
  this->filename = basename;
}

/* -------------------------------------------------------------------------- */
bool dumper::toVTKConnectivity::write_reorder_initialized = false;
std::map<ElementType, Vector<Idx>> dumper::toVTKConnectivity::write_reorder;
/* -------------------------------------------------------------------------- */
bool dumper::toAkantuConnectivity::write_reorder_initialized = false;
std::map<ElementType, Vector<Idx>> dumper::toAkantuConnectivity::write_reorder;
/* -------------------------------------------------------------------------- */
} // namespace akantu
