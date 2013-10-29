/**
 * @file   dumper_paraview.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 26 21:52:40 2012
 *
 * @brief  implementations of DumperParaview
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
#include <io_helper.hh>

#include "dumper_paraview.hh"
#include "static_communicator.hh"
#include <fstream>

__BEGIN_AKANTU__

DumperParaview::DumperParaview(const std::string & filename,
			       const std::string & directory,
			       bool parallel) : DumperIOHelper() {
  iohelper::DumperParaview * dumper_para = new iohelper::DumperParaview();
  dumper = dumper_para;
  setBaseName(filename);

  this->setParallelContext(parallel);

  dumper_para->setMode(iohelper::BASE64);
  dumper_para->setPrefix(directory);
  dumper_para->init();
}

/* -------------------------------------------------------------------------- */
DumperParaview::~DumperParaview() {
}


/* -------------------------------------------------------------------------- */
void DumperParaview::setBaseName(const std::string & basename) {
  DumperIOHelper::setBaseName(basename);
  static_cast<iohelper::DumperParaview*>(dumper)->setVTUSubDirectory(filename + "-VTU");
}

__END_AKANTU__
