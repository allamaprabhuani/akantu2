/**
 * @file   test_structural_mechanics_model_bernoulli_beam_2.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu May 12 16:34:09 2011
 *
 * @brief
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
#include "mesh_io_msh_struct.hh"
#include "structural_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#define TYPE _bernoulli_beam_2

using namespace akantu;

//Linear load function
static void lin_load(double * position, double * load,
		     __attribute__ ((unused)) Real * normal, __attribute__ ((unused)) UInt surface_id){
  memset(load,0,sizeof(Real)*3);
  load[1]= -60000000;
 }

int main(int argc, char *argv[]){

  initialize(argc, argv);
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);

  /* -------------------------------------------------------------------------- */
  // Defining the mesh

  UInt nb_nodes=101;
  UInt nb_element=nb_nodes-1;
  UInt nb_nodes_h=101;
  UInt nb_nodes_v= nb_nodes-nb_nodes_h;

  Array<Real> & nodes = const_cast<Array<Real> &>(beams.getNodes());
  nodes.resize(nb_nodes);

  beams.addConnectivityType(_bernoulli_beam_2);
  Array<UInt> & connectivity = const_cast<Array<UInt> &>(beams.getConnectivity(_bernoulli_beam_2));
  connectivity.resize(nb_element);

  for(UInt i=0; i<nb_nodes_h; ++i) {

    nodes(i,0)=2./((Real)nb_nodes_h-1)*i;
    nodes(i,1)=0;
  }
  for(UInt i=nb_nodes_h; i<nb_nodes; ++i) {

    nodes(i,0)=2;
    nodes(i,1)=2./((Real)nb_nodes_v)*(i-nb_nodes_h);
  }

  /**nodes(0,0)=0;
     nodes(1,0)=10;
     nodes(2,0)=18;
  */

  for(UInt i=0; i<nb_element; ++i) {
    connectivity(i,0)=i;
    connectivity(i,1)=i+1;
  }

  akantu::StructuralMechanicsModel model(beams);
  model.initFull();

  StructuralMaterial mat1;
  mat1.E=2.05e11;
  mat1.I=0.00128;
  mat1.A=0.01;
  model.addMaterial(mat1);

  Array<bool> & boundary = model.getBoundary();
  Array<UInt> element_material = model.getElementMaterial(_bernoulli_beam_2);
  model.computeForcesFromFunction<_bernoulli_beam_2>(lin_load, akantu::_bft_traction);

  boundary(0,0) = true;
  boundary(0,1) = true;
  boundary(nb_nodes-1,1) = true;

  Real error;
  model.assembleStiffnessMatrix();
  model.getStiffnessMatrix().saveMatrix("Kb.mtx");

  UInt count = 0;
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

  model.getStiffnessMatrix().saveMatrix("Ka.mtx");
  model.dump();
}

