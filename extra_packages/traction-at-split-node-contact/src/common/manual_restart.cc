/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "manual_restart.hh"
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

void dumpArray(const Array<Real> & array, const std::string & fname) {
  std::ofstream outFile;
  outFile.open(fname.c_str());
  outFile.precision(9);
  outFile.setf(std::ios::scientific);

  auto size = array.size();
  auto nb_component = array.getNbComponent();
  outFile << size << std::endl;
  outFile << nb_component << std::endl;

  for (auto && v : make_view(array, nb_component)) {
    for (Int c = 0; c < nb_component; ++c) {
      if (c != 0) {
        outFile << " ";
      }
      outFile << (v)(c);
    }
    outFile << std::endl;
  }
  outFile.close();
}

void loadArray(Array<Real> & array, const std::string & fname) {
  std::ifstream inFile;
  inFile.open(fname.c_str());
  inFile.precision(9);
  inFile.setf(std::ios::scientific);
  Int size(0), nb_comp(0);
  inFile >> size;
  inFile >> nb_comp;
  AKANTU_DEBUG_ASSERT(array.getNbComponent() == nb_comp,
                      "BAD NUM OF COMPONENTS");
  AKANTU_DEBUG_ASSERT(array.size() == size,
                      "loadArray: number of data points in file ("
                          << size << ") does not correspond to array size ("
                          << array.size() << ")!!");

  array.resize(size);

  for (auto && v : make_view(array, array.getNbComponent())) {
    for (auto && c : v) {
      inFile >> c;
    }
  }
  inFile.close();
}

/* -------------------------------------------------------------------------- */
void loadRestart(akantu::SolidMechanicsModel & model, const std::string & fname,
                 akantu::UInt prank) {

  const auto & mesh = model.getMesh();
  auto spatial_dimension = model.getMesh().getSpatialDimension();
  auto & dof_manager = dynamic_cast<DOFManagerDefault &>(model.getDOFManager());
  if (prank == 0) {
    akantu::Array<akantu::Real> full_reload_array(mesh.getNbGlobalNodes(),
                                                  spatial_dimension);
    loadArray(full_reload_array, fname);
    dof_manager.getSynchronizer().scatter(model.getDisplacement(),
                                          full_reload_array);
  } else {
    dof_manager.getSynchronizer().scatter(model.getDisplacement());
  }
}

/* -------------------------------------------------------------------------- */
void loadRestart(akantu::SolidMechanicsModel & model,
                 const std::string & fname) {
  loadArray(model.getDisplacement(), fname);
}

/* -------------------------------------------------------------------------- */
void dumpRestart(akantu::SolidMechanicsModel & model, const std::string & fname,
                 akantu::UInt prank) {

  const akantu::Mesh & mesh = model.getMesh();
  const akantu::Int spatial_dimension = model.getMesh().getSpatialDimension();
  auto & dof_manager = dynamic_cast<DOFManagerDefault &>(model.getDOFManager());

  if (prank == 0) {
    akantu::Array<akantu::Real> full_array(mesh.getNbGlobalNodes(),
                                           spatial_dimension);
    dof_manager.getSynchronizer().gather(model.getDisplacement(), full_array);
    dumpArray(full_array, fname);
  } else {
    dof_manager.getSynchronizer().gather(model.getDisplacement());
  }
}

/* -------------------------------------------------------------------------- */
void dumpRestart(akantu::SolidMechanicsModel & model,
                 const std::string & fname) {
  dumpArray(model.getDisplacement(), fname);
}
