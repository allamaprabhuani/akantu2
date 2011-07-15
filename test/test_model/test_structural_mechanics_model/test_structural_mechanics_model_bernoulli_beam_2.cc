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

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"

void paraviewInit(Dumper & dumper, const StructuralMechanicsModel & model);
void paraviewDump(Dumper & dumper);
#endif

//Linear load function
static void lin_load(double * position, double * load,
		     __attribute__ ((unused)) Real * normal, __attribute__ ((unused)) UInt surface_id){
  memset(load,0,sizeof(Real)*3);
  load[1]= -60000000;
 }

int main(int argc, char *argv[]){

  initialize(&argc, &argv);
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);

  /* -------------------------------------------------------------------------- */
  // Defining the mesh

  UInt nb_nodes=101;
  UInt nb_element=nb_nodes-1;
  UInt nb_nodes_h=101;
  UInt nb_nodes_v= nb_nodes-nb_nodes_h;

    Vector<Real> & nodes = const_cast<Vector<Real> &>(beams.getNodes());
    nodes.resize(nb_nodes);

    beams.addConnecticityType(_bernoulli_beam_2);
    Vector<UInt> & connectivity = const_cast<Vector<UInt> &>(beams.getConnectivity(_bernoulli_beam_2));
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
    akantu::MeshIOMSHStruct mesh_io;
    mesh_io.write("b_beam_2.msh", beams);

/* -------------------------------------------------------------------------- */
    // Defining the forces


  akantu::ElementType type = akantu::_bernoulli_beam_2;

  akantu::StructuralMechanicsModel * model;

  model = new akantu::StructuralMechanicsModel(beams);
  model->initModel();

  StructuralMaterial mat1;
  mat1.E=2.05e11;
  mat1.I=0.00128;
  mat1.A=0.01;


  model->addMaterial(mat1);

  model->initVectors();
  model->initImplicitSolver();

  const Real M = 3600; // Momentum at 3

  Vector<Real> & forces = model->getForce();
  Vector<Real> & displacement = model->getDisplacement();
  Vector<bool> & boundary = model->getBoundary();

  Vector<UInt> element_material = model->getElementMaterial(_bernoulli_beam_2);

  forces.clear();
  displacement.clear();


  model->computeForcesFromFunction(lin_load, akantu::_bft_forces);

  /**forces(0,2)=-N;
  forces(nb_nodes-1,2)=N;

  for (UInt i = 0; i < nb_nodes; ++i) {
    forces(i,1)=-N/nb_nodes;
  }
  */

  boundary(0,0) = true;
  boundary(0,1) = true;
  //boundary(0,2) = true;
  //boundary(nb_nodes-1,0) = true;
  boundary(nb_nodes-1,1) = true;
  //boundary(nb_nodes-1,2) = true;

  Real error;

  model->assembleStiffnessMatrix();
  model->getStiffnessMatrix().saveMatrix("Kb.mtx");
  UInt count = 0;

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  paraviewInit(dumper, *model);
#endif

  do {
    if(count != 0) std::cerr << count << " - " << error << std::endl;
    model->updateResidual();
    model->solve();
    count++;
  } while (!model->testConvergenceIncrement(1e-10, error) && count < 10);
  std::cerr << count << " - " << error << std::endl;

  model->getStiffnessMatrix().saveMatrix("Ka.mtx");
#ifdef AKANTU_USE_IOHELPER
  paraviewDump(dumper);
#endif
}



/* -------------------------------------------------------------------------- */
/* Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
template <ElementType type> UInt paraviewType();
template <> UInt paraviewType<_segment_2>()      { return LINE1; };
template <> UInt paraviewType<_segment_3>()      { return LINE2; };
template <> UInt paraviewType<_triangle_3>()     { return TRIANGLE1; };
template <> UInt paraviewType<_triangle_6>()     { return TRIANGLE2; };
template <> UInt paraviewType<_quadrangle_4>()   { return QUAD1; };
template <> UInt paraviewType<_tetrahedron_4>()  { return TETRA1; };
template <> UInt paraviewType<_tetrahedron_10>() { return TETRA2; };
template <> UInt paraviewType<_hexahedron_8>()   { return HEX1; };
template <> UInt paraviewType<_bernoulli_beam_2>(){ return LINE1; };
/* -------------------------------------------------------------------------- */
void paraviewInit(Dumper & dumper, const StructuralMechanicsModel & model) {
  UInt spatial_dimension = ElementClass<TYPE>::getSpatialDimension();
  UInt nb_nodes   = model.getFEM().getMesh().getNbNodes();
  UInt nb_element = model.getFEM().getMesh().getNbElement(TYPE);

  std::stringstream filename; filename << "beam";

  dumper.SetMode(TEXT);
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, filename.str().c_str());
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(TYPE).values,
			 paraviewType<TYPE>(), nb_element, C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values,
			  3, "displacements");
  dumper.AddNodeDataField(model.getForce().values,
			  3, "applied_force");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(Dumper & dumper) {
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */

#endif
