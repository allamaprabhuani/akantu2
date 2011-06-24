/**
 * @file   test_heat_transfer_model_cube3d.cc
 * @author Rui WANG<rui.wang@epfl.ch>
 * @date   Tue May 17 11:31:22 2011
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

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "heat_transfer_model.hh"
#include "pbc_synchronizer.hh"
#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

void paraviewInit(akantu::HeatTransferModel * model,Dumper & dumper);
void paraviewDump(Dumper & dumper);

akantu::UInt spatial_dimension = 2;
akantu:: ElementType type = akantu::_triangle_3;
akantu::UInt paraview_type = TRIANGLE1;

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  akantu::initialize(&argc,&argv);

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;

  mesh_io.read("square_tri3.msh", mesh);
  
  akantu::HeatTransferModel * model;
  akantu::UInt nb_nodes;
  akantu::UInt nb_element;

  model = new akantu::HeatTransferModel(mesh);
  /* -------------------------------------------------------------------------- */
  model->readMaterials("material.dat");
  model->initModel();
  model->initVectors();
  model->getHeatFlux().clear();
  model->getCapacityLumped().clear();
  model->getTemperatureGradient(type).clear();
  /* -------------------------------------------------------------------------- */
  model->initPBC(1,1,1);
  model->assembleCapacityLumped(type);
  /* -------------------------------------------------------------------------- */
  nb_nodes = model->getFEM().getMesh().getNbNodes();
  nb_element = model->getFEM().getMesh().getNbElement(type);
  nb_nodes = model->getFEM().getMesh().getNbNodes();
  /* ------------------------------------------------------------------------ */
  //get stable time step
  akantu::Real time_step = model->getStableTimeStep()*0.8;
  cout<<"time step is:"<<time_step<<endl;
  model->setTimeStep(time_step);
  /* -------------------------------------------------------------------------- */
  /// boundary conditions
  const akantu::Vector<akantu::Real> & nodes = model->getFEM().getMesh().getNodes();
  akantu::Vector<bool> & boundary = model->getBoundary();
  akantu::Vector<akantu::Real> & temperature = model->getTemperature();
  akantu::Vector<akantu::Real> & heat_flux = model->getHeatFlux();
  akantu::Real eps = 1e-15;

  double t1, t2, length;
  t1 = 300.;
  t2 = 100.;
  length = 1.;

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    temperature(i) = 100.;

    akantu::Real dx = nodes(i,0) - length/4.;
    akantu::Real dy = 0.0;
    akantu::Real dz = 0.0;

    if (spatial_dimension > 1) dy = nodes(i,1) - length/4.;
    if (spatial_dimension == 3) dz = nodes(i,2) - length/4.;
    akantu::Real d = sqrt(dx*dx + dy*dy + dz*dz);
    //    if(dx < 0.0){
    if(d < 0.1){
      boundary(i) = true;
      temperature(i) = 300.;
    }
  }
  /* -------------------------------------------------------------------------- */
  DumperParaview dumper;
  paraviewInit(model,dumper);
  /* ------------------------------------------------------------------------ */
  // //for testing
  int max_steps = 100000;
  /* ------------------------------------------------------------------------ */
  for(int i=0; i<max_steps; i++)
    {
      model->updateHeatFlux();
      model->updateTemperature();
     
      if(i % 100 == 0)
	paraviewDump(dumper);
      if(i % 10 == 0)
      std::cout << "Step " << i << "/" << max_steps << std::endl;
    }
  cout<< "\n\n Stable Time Step is : " << time_step << "\n \n" <<endl;

  return 0;
}
/* -------------------------------------------------------------------------- */

void paraviewInit(akantu::HeatTransferModel * model, Dumper & dumper) {
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::UInt nb_element = model->getFEM().getMesh().getNbElement(type);
  

  dumper.SetMode(TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "coordinates2");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 paraview_type, nb_element, C_MODE);
   dumper.AddNodeDataField(model->getTemperature().values,
    1, "temperature");
  dumper.AddNodeDataField(model->getHeatFlux().values,
   			  1, "heat_flux");
  dumper.AddNodeDataField(model->getCapacityLumped().values,
   			  1, "capacity_lumped");
  dumper.AddElemDataField(model->getTemperatureGradient(type).values,
    			  spatial_dimension, "temperature_gradient");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */

void paraviewDump(Dumper & dumper) {
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
