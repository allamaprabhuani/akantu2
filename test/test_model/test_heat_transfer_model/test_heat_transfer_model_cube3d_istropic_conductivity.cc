/**
 * @file   test_heat_transfer_model_cube3d_istropic_conductivity.cc
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "heat_transfer_model.hh"
//#include "material.hh"
// basic file operations
#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

 akantu::UInt spatial_dimension = 3;
 akantu:: ElementType type = akantu::_tetrahedron_4;
 akantu::UInt paraview_type = TETRA1;



//just for checking
 


void paraviewInit(Dumper & dumper);
void paraviewDump(Dumper & dumper);

akantu::HeatTransferModel * model;

akantu::UInt nb_nodes;
akantu::UInt nb_element;


akantu::Real density;
akantu::Real conductivity[3][3];
akantu::Real capacity;


int readMaterial () {
 
  string str;
  ifstream myfile;
 

  myfile.open("Material.dat");
  if(!myfile) //Always test the file open.
  {
    cout<<"Error opening output file"<<endl;
    return -1;
  }
  
  getline(myfile, str);
  density=atof(str.c_str());
    
  getline(myfile, str);
  capacity=atof(str.c_str());
  
  getline(myfile, str);
  char * cstr, *p;
  char * tmp_cstr;
  cstr = new char [str.size()+1];
  strcpy (cstr, str.c_str());
  tmp_cstr = cstr;
  for(int i=0;i<3;i++)
    for(int j=0; j<3;j++)
      {
	p=strtok(tmp_cstr, " "); tmp_cstr = NULL;
	conductivity[i][j]= atof(p);
	cout<<conductivity[i][j]<<endl;
      }
 

   return 0;
}


int main(int argc, char *argv[])
{
  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;

  mesh_io.read("cube1.msh", mesh);
  readMaterial();

  cout<<"The density of the material is:"<< density <<endl;
  cout<<"The capacity of the material is:"<< capacity <<endl;

  model = new akantu::HeatTransferModel(mesh);
  //model initialization
  model->initModel();
  //initialize the vectors
  model->initVectors();

  nb_nodes = model->getFEM().getMesh().getNbNodes();
  nb_element = model->getFEM().getMesh().getNbElement(type);

 
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  model->getHeatFlux().clear();
  model->getLumped().clear();
  model->getTemperatureGradient(type).clear();

  // akantu::debug::setDebugLevel(akantu::dblDump);
  // std::cout << model->getTemperatureGradient(type) << std::endl;
  // akantu::debug::setDebugLevel(akantu::dblWarning);

 

  model->setDensity(density);
  model->setCapacity(capacity);
  model->SetConductivityMatrix(conductivity);
 
  //get stable time step
  akantu::Real time_step = model->getStableTimeStep()*0.8;
  
  cout<<"time step is:"<<time_step<<endl;
  model->setTimeStep(time_step);

  /// boundary conditions
  const akantu::Vector<akantu::Real> & nodes = model->getFEM().getMesh().getNodes();
  akantu::Vector<bool> & boundary = model->getBoundary();
  akantu::Vector<akantu::Real> & temperature = model->getTemperature();
  akantu::Vector<akantu::Real> & heat_flux = model->getHeatFlux();
  akantu::Real eps = 1e-15;
 
  double t1, t2, length;


  length = 1.;

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

  DumperParaview dumper;
  paraviewInit(dumper);
  model->assembleMassLumped(type);


  // //for testing
  int max_steps = 100000;

  for(int i=0; i<max_steps; i++)
    {
     
      model->updateHeatFlux();
      model->updateTemperature();
     
      if(i % 100 == 0)
	paraviewDump(dumper);
      if(i % 10000 == 0)
      std::cout << "Step " << i << "/" << max_steps << std::endl;
    }
  cout<< "\n\n Stable Time Step is : " << time_step << "\n \n" <<endl;

  return 0;
}

void paraviewInit(Dumper & dumper) {
  dumper.SetMode(TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "coordinates_cube3d_istropic_conductivity");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 paraview_type, nb_element, C_MODE);
   dumper.AddNodeDataField(model->getTemperature().values,
    1, "temperature");
  dumper.AddNodeDataField(model->getHeatFlux().values,
   			  1, "heat_flux");
  dumper.AddNodeDataField(model->getLumped().values,
   			  1, "lumped");
  dumper.AddElemDataField(model->getTemperatureGradient(type).values,
    			  spatial_dimension, "temperature_gradient");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}


void paraviewDump(Dumper & dumper) {
  dumper.Dump();
}

