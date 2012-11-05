/**
 * @file   test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date   Fri Jul 15 19:41:58 2011
 *
 * @brief  Computation of the analytical exemple 1.1 in the TGC vol 6
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

#define TYPE _bernoulli_beam_2

using namespace akantu;

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"


static void paraviewInit(iohelper::Dumper & dumper, const StructuralMechanicsModel & model);
static void paraviewDump(iohelper::Dumper & dumper);
#endif

//Linear load function
static void lin_load(double * position, double * load,
		     __attribute__ ((unused)) Real * normal, __attribute__ ((unused)) UInt surface_id){
  memset(load,0,sizeof(Real)*3);
  if (position[0]<=10){
  load[1]= -6000;
 }
}

int main(int argc, char *argv[]){

  initialize(argc, argv);
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);

  /* -------------------------------------------------------------------------- */
  // Defining the mesh

  UInt nb_nodes=401;
  UInt nb_nodes_1=200;
  UInt nb_nodes_2=nb_nodes-nb_nodes_1 - 1;
  UInt nb_element=nb_nodes-1;

  Vector<Real> & nodes = const_cast<Vector<Real> &>(beams.getNodes());
  nodes.resize(nb_nodes);

  beams.addConnectivityType(_bernoulli_beam_2);
  Vector<UInt> & connectivity = const_cast<Vector<UInt> &>(beams.getConnectivity(_bernoulli_beam_2));
  connectivity.resize(nb_element);

  for(UInt i=0; i<nb_nodes; ++i) {
     nodes(i,1)=0;
  }
   for (UInt i = 0; i < nb_nodes_1; ++i) {
     nodes(i,0)=10.*i/((Real)nb_nodes_1);

   }
   nodes(nb_nodes_1,0)=10;

   for (UInt i = 0; i < nb_nodes_2; ++i) {
     nodes(nb_nodes_1 + i + 1,0)=10+ 8.*(i+1)/((Real)nb_nodes_2);

   }

  for(UInt i=0; i<nb_element; ++i) {

    connectivity(i,0)=i;
    connectivity(i,1)=i+1;
  }
  akantu::MeshIOMSH mesh_io;
  mesh_io.write("b_beam_2.msh", beams);

  /* -------------------------------------------------------------------------- */
  // Defining the materials


  //  akantu::ElementType type = akantu::_bernoulli_beam_2;

  akantu::StructuralMechanicsModel * model;

  model = new akantu::StructuralMechanicsModel(beams);

  StructuralMaterial mat1;
  mat1.E=3e10;
  mat1.I=0.0025;
  mat1.A=0.01;


  model->addMaterial(mat1);

  StructuralMaterial mat2 ;
  mat2.E=3e10;
  mat2.I=0.00128;
  mat2.A=0.01;

  model->addMaterial(mat2);

  /* -------------------------------------------------------------------------- */
  // Defining the forces

  model->initModel();

  model->initVectors();

  const Real M = -3600; // Momentum at 3

  Vector<Real> & forces = model->getForce();
  Vector<Real> & displacement = model->getDisplacement();
  Vector<bool> & boundary = model->getBoundary();
  const Vector<Real> & N_M  = model->getStress(_bernoulli_beam_2);

  Vector<UInt> & element_material = model->getElementMaterial(_bernoulli_beam_2);

  forces.clear();
  displacement.clear();

  for (UInt i = 0; i < nb_nodes_2; ++i) {
    element_material(i+nb_nodes_1)=1;
}


  forces(nb_nodes-1,2) += M;

  model->computeForcesFromFunction(lin_load, akantu::_bft_traction);

  /* -------------------------------------------------------------------------- */
  // Defining the boundary conditions

  boundary(0,0) = true;
  boundary(0,1) = true;
  boundary(0,2) = true;
  boundary(nb_nodes_1,1) = true;
  boundary(nb_nodes-1,1) = true;
  /* -------------------------------------------------------------------------- */
  // Solve

  model->initImplicitSolver();

  Real error;

  model->assembleStiffnessMatrix();
  model->getStiffnessMatrix().saveMatrix("Kb.mtx");
  UInt count = 0;

#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, *model);
#endif

  do {
    if(count != 0) std::cerr << count << " - " << error << std::endl;
    model->updateResidual();
    model->solve();
    count++;
  } while (!model->testConvergenceIncrement(1e-10, error) && count < 10);
  std::cerr << count << " - " << error << std::endl;

  /* -------------------------------------------------------------------------- */
  // Post-Processing

  model->computeStressOnQuad();

  model->getStiffnessMatrix().saveMatrix("Ka.mtx");
  std::cout<< " d1 = " << displacement(nb_nodes_1,2) << std::endl;
  std::cout<< " d2 = " << displacement(nb_nodes-1,2) << std::endl;
  std::cout<< " M1 = " << N_M(0,1) << std::endl;
  std::cout<< " M2 = " << N_M(2*(nb_nodes-2),1) << std::endl;


#ifdef AKANTU_USE_IOHELPER
  paraviewDump(dumper);
#endif
}



/* -------------------------------------------------------------------------- */
/* iohelper::Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
// template <ElementType type> static iohelper::ElemType paraviewType();
// template <> iohelper::ElemType paraviewType<_segment_2>()      { return iohelper::LINE1; }
// template <> iohelper::ElemType paraviewType<_segment_3>()      { return iohelper::LINE2; }
// template <> iohelper::ElemType paraviewType<_triangle_3>()     { return iohelper::TRIANGLE1; }
// template <> iohelper::ElemType paraviewType<_triangle_6>()     { return iohelper::TRIANGLE2; }
// template <> iohelper::ElemType paraviewType<_quadrangle_4>()   { return iohelper::QUAD1; }
// template <> iohelper::ElemType paraviewType<_tetrahedron_4>()  { return iohelper::TETRA1; }
// template <> iohelper::ElemType paraviewType<_tetrahedron_10>() { return iohelper::TETRA2; }
// template <> iohelper::ElemType paraviewType<_hexahedron_8>()   { return iohelper::HEX1; }
// template <> iohelper::ElemType paraviewType<_bernoulli_beam_2>(){ return iohelper::BEAM2; }
/* -------------------------------------------------------------------------- */
void paraviewInit(iohelper::Dumper & dumper, const StructuralMechanicsModel & model) {
  // UInt spatial_dimension = ElementClass<TYPE>::getSpatialDimension();
  // UInt nb_nodes   = model.getFEM().getMesh().getNbNodes();
  // UInt nb_element = model.getFEM().getMesh().getNbElement(TYPE);

  // std::stringstream filename; filename << "beam";

#pragma message "To change with new dumper"
  // dumper.SetMode(iohelper::TEXT);
  // dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
  // 		   spatial_dimension, nb_nodes, filename.str().c_str());
  // dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(TYPE).values,
  // 			 paraviewType<TYPE>(), nb_element, iohelper::C_MODE);
  // dumper.AddNodeDataField(model.getDisplacement().values,
  // 			  3, "displacements");
  // dumper.AddNodeDataField(model.getForce().values,
  // 			  3, "applied_force");
  // dumper.AddElemDataField(model.getStress(_bernoulli_beam_2).values,
  // 			  2, "stress");
  // dumper.SetPrefix("paraview/");
  // dumper.Init();
  // dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(iohelper::Dumper & dumper) {
  //  dumper.Dump();
}

/* -------------------------------------------------------------------------- */

#endif
