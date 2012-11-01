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
#include "dumper_paraview.hh"
#include "static_communicator.hh"

__BEGIN_AKANTU__

DumperParaview::DumperParaview(const std::string & filename,
			       const std::string & directory) : DumperIOHelper() {
  this->filename = filename;
  iohelper::DumperParaview * dumper_para = new iohelper::DumperParaview();
  dumper = dumper_para;

  UInt whoami = StaticCommunicator::getStaticCommunicator().whoAmI();
  UInt nproc  = StaticCommunicator::getStaticCommunicator().getNbProc();
  dumper_para->setMode(iohelper::TEXT);
  dumper_para->setParallelContext(whoami, nproc);

  dumper_para->setVTUSubDirectory(filename + "-VTU");
  dumper_para->setPrefix(directory);
  dumper_para->init();
}

__END_AKANTU__
