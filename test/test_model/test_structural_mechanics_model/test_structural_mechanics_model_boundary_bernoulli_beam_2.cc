/**
 * @file   test_structural_mechanics_model_boundary_bernoulli_beam_2.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu May 26 12:52:54 2011
 *
 * @brief  Test the computation of linear load
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
/* -------------------------------------------------------------------------- */

#define TYPE _bernoulli_beam_2

using namespace akantu;

static void lin_load(double * position, double * load,
		     __attribute__ ((unused)) Real * normal, __attribute__ ((unused)) UInt surface_id){
  memset(load,0,sizeof(Real)*3);
  //load[1]= (position[0])*200000000;
  load[1]= 200000000;
  }

int main(int argc, char *argv[]){

  initialize(argc, argv);
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);

/* -------------------------------------------------------------------------- */
  // Defining the mesh

  UInt nb_nodes=2;
  UInt nb_element=nb_nodes-1;
  UInt nb_nodes_h=2;
  UInt nb_nodes_v= nb_nodes-nb_nodes_h;

    Vector<Real> & nodes = const_cast<Vector<Real> &>(beams.getNodes());
    nodes.resize(nb_nodes);

    beams.addConnecticityType(_bernoulli_beam_2);
    Vector<UInt> & connectivity = const_cast<Vector<UInt> &>(beams.getConnectivity(_bernoulli_beam_2));
    connectivity.resize(nb_element);

    for(UInt i=0; i<nb_nodes_h; ++i) {

      nodes(i,0)=2./((Real)nb_nodes_h)*i;
      nodes(i,1)=0;
    }
    for(UInt i=nb_nodes_h; i<nb_nodes; ++i) {

      nodes(i,0)=2;
      nodes(i,1)=2./((Real)nb_nodes_v)*(i-nb_nodes_h);
    }

    for(UInt i=0; i<nb_element; ++i) {

      connectivity(i,0)=i;
      connectivity(i,1)=i+1;
    }
    akantu::MeshIOMSH mesh_io;
    mesh_io.write("b_beam_2.msh", beams);

/* -------------------------------------------------------------------------- */
    // Defining the forces


    //  akantu::ElementType type = akantu::_bernoulli_beam_2;

  akantu::StructuralMechanicsModel * model;

  model = new akantu::StructuralMechanicsModel(beams);

  model->initModel();
  model->initVectors();

  Vector<Real> & forces = model->getForce();

  forces.clear();

  model->computeForcesFromFunction(lin_load, akantu::_bft_traction);


}
