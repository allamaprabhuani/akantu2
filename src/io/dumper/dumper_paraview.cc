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
			       const std::string & directory,
                               bool parallel) : DumperIOHelper() {

  iohelper::DumperParaview * dumper_para = new iohelper::DumperParaview();
  dumper = dumper_para;
  setBaseName(filename);

  this->setParallelContext(parallel);

  dumper_para->setMode(iohelper::BASE64);
  dumper_para->setPrefix(directory);
  dumper_para->init();

  //  this->directory = dumper_para->getPrefix();

  // current_time = 0.;
  // time_step = 0.;
}

DumperParaview::~DumperParaview() {
  // if(pvd_file.is_open()) {
  //   pvd_file << "  </Collection>" << std::endl
  // 	     << "</VTKFile>" << std::endl;
  //   pvd_file.close();
  // }
}


// void DumperParaview::dump() {
//   DumperIOHelper::dump();

//   if(!pvd_file.is_open()) {
//     pvd_file.open(directory + filename + ".pvd");

//     if(!pvd_file.good()) {
//       AKANTU_EXCEPTION("DumperParaview was not able to open the file \"" << directory << filename << ".pvd\"");
//     }

//     pvd_file << "<?xml version=\"1.0\"?>" << std::endl
// 	     << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
// 	     << "  <Collection>" << std::endl;
//   }

//   if(pvd_file.is_open()) {
//     pvd_file << "    <DataSet timestep=\"" << current_time << "\" group=\"\" part=\"0\" file=\""
// 	     << last_filename << ".pvtu\"/>" << std::endl;
//   }
//   current_time += time_step;
// }

__END_AKANTU__
 
