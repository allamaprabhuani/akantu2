/**
 * @file   test_heat_transfer_model_cube3d_istropic_conductivity.cc
 *
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date   Sun May 01 19:14:43 2011
 *
 * @brief  test of the class HeatTransferModel on the 3d cube
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "heat_transfer_model.hh"
#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;

/* -------------------------------------------------------------------------- */
akantu::UInt spatial_dimension = 3;
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;

  mesh_io.read("cube1.msh", mesh);

  akantu::HeatTransferModel model(mesh);
  //initialize everything
  model.initFull("material.dat");
  //assemble the lumped capacity
  model.assembleCapacityLumped();

  //get stable time step
  akantu::Real time_step = model.getStableTimeStep()*0.8;
  cout<<"time step is:"<<time_step<<endl;
  model.setTimeStep(time_step);

  /// boundary conditions
  const akantu::Vector<akantu::Real> & nodes = model.getFEM().getMesh().getNodes();
  akantu::Vector<bool> & boundary = model.getBoundary();
  akantu::Vector<akantu::Real> & temperature = model.getTemperature();
  akantu::Real eps = 1e-15;

  double length = 1.;
  akantu::UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    //temperature(i) = t1 - (t1 - t2) * sin(nodes(i, 0) * M_PI / length);
    temperature(i) = 100.;

    if(nodes(i,0) < eps) {
      boundary(i) = true;
      temperature(i) = 300.;
    }
    //set the second boundary condition
    if(std::abs(nodes(i,0) - length) < eps) {
      boundary(i) = true;
      temperature(i) = 300.;
    }
    // //to insert a heat source
    //  if(std::abs(nodes(i,0) - length/2.) < 0.025 && std::abs(nodes(i,1) - length/2.) < 0.025 && std::abs(nodes(i,2) - length/2.) < 0.025) {
    //   boundary(i) = true;
    //  temperature(i) = 300.;
    //  }
  }

  model.updateResidual();
  model.setBaseName("heat_transfer_cube3d_istropic_conductivity");
  model.addDumpField("temperature"     );
  model.addDumpField("temperature_rate");
  model.addDumpField("residual"        );
  model.addDumpField("capacity_lumped" );
  model.dump();

  // //for testing
  int max_steps = 1000;

  for(int i=0; i<max_steps; i++)
    {
      model.explicitPred();
      model.updateResidual();
      model.explicitCorr();

      if(i % 100 == 0) model.dump();
      if(i % 10000 == 0)
      std::cout << "Step " << i << "/" << max_steps << std::endl;
    }
  cout<< "\n\n Stable Time Step is : " << time_step << "\n \n" <<endl;

  return 0;
}
