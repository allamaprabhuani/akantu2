/**
 * @file   test_structural_mechanics_model_bernoulli_beam_local_force.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Fri Jan  6 12:31:43 2012
 *
 * @brief  Test force defined locally
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
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "structural_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#define TYPE _bernoulli_beam_3

using namespace akantu;

//Linear load function
static void lin_load(double * position, double * load,
		     __attribute__ ((unused)) Real * normal,
		     __attribute__ ((unused)) UInt surface_id){
  memset(load,0,sizeof(Real)*6);
  load[1]= 60.;
  }

int main(int argc, char *argv[]){
  initialize(argc, argv);
  Mesh beams(3);
  debug::setDebugLevel(dblWarning);
  std::cout<<"Initialisation"<<std::endl;
  /* -------------------------------------------------------------------------- */
  // Defining the mesh
  
  UInt nb_nodes=3;
  UInt nb_element=nb_nodes-1;

  Vector<Real> & nodes = const_cast<Vector<Real> &>(beams.getNodes());
  nodes.resize(nb_nodes);

  beams.addConnectivityType(_bernoulli_beam_3);
  Vector<UInt> & connectivity = const_cast<Vector<UInt> &>(beams.getConnectivity(_bernoulli_beam_3));

  connectivity.resize(nb_element);

  beams.initNormals();

  Vector<Real> & normals = const_cast<Vector<Real> &>(beams.getNormals(_bernoulli_beam_3));
  normals.resize(nb_element);

  nodes(0,0)=0.;
  nodes(0,1)=0.;
  nodes(0,2)=0.;

  nodes(1,0)=1.;
  nodes(1,1)=0.;
  nodes(1,2)=0.;

  nodes(2,0)=1.;
  nodes(2,1)=1.;
  nodes(2,2)=0.;
  
  for (UInt i=0; i<nb_element; ++i) {
    
    connectivity(i,0)=i;
    connectivity(i,1)=i+1;

    normals(i,0)=0.;
    normals(i,1)=0.;
    normals(i,2)=1.;
  }

  normals(1,0)=1.;
  normals(1,1)=0.;
  normals(1,2)=0.;

  akantu::MeshIOMSH mesh_io;
  mesh_io.write("beam_local_force.msh", beams);
  std::cout<<"Mesh definition"<<std::endl;
  /* -------------------------------------------------------------------------- */
  // Defining the materials


  akantu::StructuralMechanicsModel model(beams);

  StructuralMaterial mat1;
  mat1.E=2.05e11;
  mat1.Iz=2.8e-9;
  mat1.A=1.33e-4;
  mat1.Iy=2.8e-9;
  mat1.GJ=220.77;

  model.addMaterial(mat1);

  std::cout<<"Material Definition"<<std::endl;
  /* -------------------------------------------------------------------------- */
  // Defining the forces

  model.initFull();

  Vector<bool> & boundary = model.getBoundary();

  model.computeForcesFromFunction<_bernoulli_beam_3>(lin_load, akantu::_bft_stress);
  std::cout<<"Force Definition"<<std::endl;
 /* -------------------------------------------------------------------------- */
  // Defining the boundary conditions

  boundary(0,0) = true;
  boundary(0,1) = true;
  boundary(0,2) = true;
  boundary(0,3) = true;
  boundary(0,4) = true;
  boundary(0,5) = true;
  std::cout<<"BC Definition"<<std::endl;
  /* -------------------------------------------------------------------------- */
  // Solve
  Real error;

  model.assembleStiffnessMatrix();
  std::cout<<"Assemble Done"<<std::endl;
  model.getStiffnessMatrix().saveMatrix("Klocb.mtx");
  UInt count = 0;
  std::cout<<"Matrix saved"<<std::endl;

  model.addDumpField("displacememt");
  model.addDumpField("rotation");
  model.addDumpField("force");
  model.addDumpField("momentum");

  do {
    if(count != 0) std::cerr << count << " - " << error << std::endl;
    model.updateResidual();
    model.solve();
    count++;
  } while (!model.testConvergenceIncrement(1e-10, error) && count < 10);
  std::cerr << count << " - " << error << std::endl;

  /* -------------------------------------------------------------------------- */
  // Post-Processing

  model.computeStresses();

  model.getStiffnessMatrix().saveMatrix("Kloca.mtx");
 
  model.dump();
}
