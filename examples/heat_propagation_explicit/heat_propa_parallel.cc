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
#include "mesh_partition_scotch.hh"
#ifdef AKANTU_USE_QVIEW
#include <libqview.h>
#endif //AKANTU_USE_QVIEW
/* -------------------------------------------------------------------------- */
using namespace std;
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

void paraviewInit(const string & name,akantu::HeatTransferModel & model,Dumper & dumper,
		  const akantu::GhostType & ghost_type);
void ecinDump(akantu::HeatTransferModel & model,akantu::UInt current_step);

akantu::UInt spatial_dimension = 3;
akantu:: ElementType type = akantu::_tetrahedron_4;
akantu::UInt paraview_type = TETRA1;

std::string case_name;
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  akantu::initialize(&argc,&argv);

  case_name = argv[4];
  akantu::Mesh mesh(spatial_dimension);
  akantu::StaticCommunicator * comm =
    akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::debug::setDebugLevel(akantu::dblWarning);
  std::stringstream filename;
  filename << "log" << prank;
  akantu::debug::setLogFile(filename.str());
  akantu::MeshPartition * partition = NULL;

  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read(argv[3], mesh);
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }
  akantu::HeatTransferModel model(mesh);
  model.initParallel(partition);


  mesh.computeBoundingBox();

  /* -------------------------------------------------------------------------- */
  model.readMaterials("material.dat");
  model.initModel();
  model.initVectors();

  model.getTemperature().clear();
  model.getTemperatureRate().clear();
  /* -------------------------------------------------------------------------- */
  model.assembleCapacityLumped();
  /* -------------------------------------------------------------------------- */
  akantu::UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  /* ------------------------------------------------------------------------ */
  //get stable time step
  akantu::Real time_step = model.getStableTimeStep()*0.8;
  if (prank == 0)
    std::cerr << "time step is:"<<time_step << std::endl;
  //  model.setTimeStep(time_step);
  if (atof(argv[2]) > time_step)
    AKANTU_DEBUG_ERROR("the required timestep is higher than the stable timestep of the mesh");
  model.setTimeStep(atof(argv[2]));

  /* -------------------------------------------------------------------------- */
  /// boundary conditions
  const akantu::Vector<akantu::Real> & nodes = model.getFEM().getMesh().getNodes();
  akantu::Vector<bool> & boundary = model.getBoundary();
  akantu::Vector<akantu::Real> & temperature = model.getTemperature();

  double t1, t2;
  t1 = 150.;
  t2 = 100.;

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if (mesh.isPureGhostNode(i)) continue;

    temperature(i) = t2;
    akantu::Real dz = nodes(i,2) - mesh.getZMin();
    //    akantu::Real size = mesh.getZMax() - mesh.getZMin();
    if(fabs(dz) < 0.1){
      boundary(i) = true;
      temperature(i) = t1;
    }
    // else {
    //   temperature(i) = t1-(t1-t2)*dz/size;
    // }
  }

  /* -------------------------------------------------------------------------- */
  DumperParaview dumper;
  DumperParaview dumper_ghost;

  model.updateResidual();
  
  paraviewInit(case_name,model,dumper,akantu::_not_ghost);
  //  paraviewInit("heat-test",model,dumper_ghost,akantu::_ghost);

  /* -------------------------------------------------------------------------- */
  int max_steps = atoi(argv[1]);
  /* ------------------------------------------------------------------------ */
#ifdef AKANTU_USE_QVIEW
  QView qv;
  qv.setMode(NET_MODE);
  qv.initLibQview(prank);
  qv.beginTask("main",max_steps);
#endif

  /* ------------------------------------------------------------------------ */
  for(int i=0; i<max_steps; i++)
    {
#ifdef AKANTU_USE_QVIEW
      qv.setCurrentStep(i);
#else
      if(i % 100 == 0){
       	if(prank == 0)
      	  std::cout << "Step " << i << "/" << max_steps << std::endl;
      }
#endif

      model.explicitPred();
      model.updateResidual();
      model.solveExplicitLumped();
      model.explicitCorr();

#ifdef AKANTU_USE_IOHELPER
      if(i % 1000 == 0){
	dumper.Dump();
	//	dumper_ghost.Dump();
      }
      ecinDump(model,i);
#endif
    }
#ifdef AKANTU_USE_QVIEW
  qv.endTask();
#endif
  akantu::finalize();
  return 0;
}

/* -------------------------------------------------------------------------- */
void paraviewInit(const string & name, akantu::HeatTransferModel & model, Dumper & dumper,
		       const akantu::GhostType & ghost_type) {
  akantu::UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  akantu::UInt nb_element = model.getFEM().getMesh().getNbElement(type,ghost_type);

  akantu::StaticCommunicator * comm =
    akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  stringstream str;
  str << name;
  if (ghost_type == akantu::_ghost)
    str << "_ghost";

  //  dumper.SetMode(TEXT);
  dumper.SetParallelContext(prank, psize);
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, str.str().c_str());

  if (nb_element)
    dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(type,ghost_type).values,
			   paraview_type, nb_element, C_MODE);
  else
    dumper.SetConnectivity(NULL,paraview_type, nb_element, C_MODE);
  dumper.AddNodeDataField(model.getTemperature().values,
			  1, "temperature");
  // dumper.AddNodeDataField(model.getResidual().values,
  //  			  1, "residual");
  // dumper.AddNodeDataField(model.getCapacityLumped().values,
  //  			  1, "capacity_lumped");
  if (nb_element){
    dumper.AddElemDataField(model.getTemperatureGradient(type, ghost_type).values,
			      spatial_dimension, "temperature_gradient");

    // dumper.AddElemDataField(model.getTemperatureOnQpoints(type, ghost_type).values,
    // 			      1, "temperature_qpoints");

    dumper.AddElemDataField(model.getConductivityOnQpoints(type, ghost_type).values,
			      spatial_dimension*spatial_dimension, "conductivity");
  }
  dumper.SetPrefix("./");
  dumper.Init();
}
/* -------------------------------------------------------------------------- */

void ecinDump(akantu::HeatTransferModel & model,akantu::UInt current_step) {
  
  const akantu::Mesh & mesh = model.getFEM().getMesh();
  akantu::UInt nb_nodes = mesh.getNbNodes();
  const akantu::Vector<akantu::Real> & nodes = mesh.getNodes();
  akantu::Vector<bool> & boundary = model.getBoundary();
  akantu::Vector<akantu::Real> & temperature = model.getTemperature();

  akantu::Real temp_avg = 0;
  akantu::UInt nb = 0;

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if (mesh.isPureGhostNode(i)) continue;

    akantu::Real dz = mesh.getZMax() - nodes(i,2);
    if(fabs(dz) < 0.1){
      temp_avg += temperature(i);
      ++nb;
    }
  }

  akantu::StaticCommunicator * comm = 
    akantu::StaticCommunicator::getStaticCommunicator();
  
  comm->allReduce(&temp_avg,1,akantu::_so_sum);
  comm->allReduce(&nb,1,akantu::_so_sum);

  temp_avg /= nb;
  
  akantu::Int prank = comm->whoAmI();
  if (prank == 0){
    std::ofstream file_out;
    std::stringstream fname;
    fname << "ecin-" << case_name << ".plot";

    if (current_step == 0)
      file_out.open(fname.str().c_str(),ios_base::trunc);
    else
      file_out.open(fname.str().c_str(),ios_base::app);

    file_out << current_step*model.getTimeStep() << " " << temp_avg << " " << nb << std::endl;
  }
}
