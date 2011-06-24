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

akantu::UInt spatial_dimension = 3;
akantu:: ElementType type = akantu::_tetrahedron_4;
akantu::UInt paraview_type = TETRA1;

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  akantu::initialize(&argc,&argv);

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;

  mesh_io.read("double_cube_tet4.msh", mesh);
  
  akantu::HeatTransferModel * model;
  akantu::UInt nb_nodes;
  akantu::UInt nb_element;

  mesh.computeBoundingBox();
  std::map<akantu::UInt,akantu::UInt> pbc_pair;
  akantu::MeshUtils::computePBCMap(mesh,0,pbc_pair);
  akantu::MeshUtils::computePBCMap(mesh,1,pbc_pair);

  {
    std::map<akantu::UInt,akantu::UInt>::iterator it = pbc_pair.begin();
    std::map<akantu::UInt,akantu::UInt>::iterator end = pbc_pair.end();
    
    akantu::Real * coords = mesh.getNodes().values;
    akantu::UInt dim = mesh.getSpatialDimension();
    while(it != end){
      akantu::UInt i1 = (*it).first;
      akantu::UInt i2 = (*it).second;
      
      AKANTU_DEBUG_INFO("pairing " << i1 << "(" 
			<< coords[dim*i1] << "," << coords[dim*i1+1] << "," 
			<< coords[dim*i1+2]
			<< ") with"
			<< i2 << "(" 
			<< coords[dim*i2] << "," << coords[dim*i2+1] << "," 
			<< coords[dim*i2+2]
			<< ")");	
      ++it;
    }
  }

  akantu::PBCSynchronizer synch(pbc_pair);

  model = new akantu::HeatTransferModel(mesh);
  model->createSynchronizerRegistry(model);
  model->getSynchronizerRegistry().registerSynchronizer(synch,akantu::_gst_htm_capacity);
  model->getSynchronizerRegistry().registerSynchronizer(synch,akantu::_gst_htm_temperature);

  /* -------------------------------------------------------------------------- */
  model->readMaterials("material.dat");
  model->initModel();
  model->initVectors();
  model->getHeatFlux().clear();
  model->getCapacityLumped().clear();
  model->getTemperatureGradient(type).clear();
  /* -------------------------------------------------------------------------- */
  model->changeLocalEquationNumberforPBC(pbc_pair,1);
  model->assembleCapacityLumped(type);
  model->getSynchronizerRegistry().synchronize(akantu::_gst_htm_capacity);
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
    akantu::Real dz = nodes(i,2) - mesh.getZMin();
    if(fabs(dz) < 0.1){
      boundary(i) = true;
      temperature(i) = 150.;
    }
  }
  /* -------------------------------------------------------------------------- */
  DumperParaview dumper;
  paraviewInit(model,dumper);
  /* ------------------------------------------------------------------------ */
  // //for testing
  int max_steps = 1000000;
  /* ------------------------------------------------------------------------ */
  for(int i=0; i<max_steps; i++)
    {
     
      model->updateHeatFlux();
      model->updateTemperature();
      model->getSynchronizerRegistry().synchronize(akantu::_gst_htm_temperature);
     
      if(i % 1000 == 0){
	paraviewDump(dumper);
	std::cout << "Step " << i << "/" << max_steps << std::endl;
      }
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
  // dumper.AddElemDataField(model->getTemperatureGradient(type).values,
  //   			  spatial_dimension, "temperature_gradient");
  // dumper.AddElemDataField(model->getbtkgt().values,
  //   			  4, "btkgt");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */

void paraviewDump(Dumper & dumper) {
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
