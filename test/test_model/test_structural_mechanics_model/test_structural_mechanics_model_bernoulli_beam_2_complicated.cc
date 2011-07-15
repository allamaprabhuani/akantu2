/**
 * @file   test_structural_mechanics_model_bernoulli_beam_2_complicated.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Jun  1 16:06:45 2011
 *
 * @brief  A very complicated structure
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
  if(position[1]>=0.-Math::tolerance) {
    if ((position[0]<=10.)){
      load[1]= -100;
    } else if (position[0]<=20.){
      load[1]= -70;
    }
  }
}

int main(int argc, char *argv[]){

  initialize(&argc, &argv);
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);

  /* -------------------------------------------------------------------------- */
  // Defining the mesh

  akantu::MeshIOMSHStruct mesh_io;
  mesh_io.read("complicated.msh", beams);
  mesh_io.write("complicated_tata.msh", beams);

  /* -------------------------------------------------------------------------- */
  // Defining the material

  const akantu::ElementType type = akantu::_bernoulli_beam_2;

  akantu::StructuralMechanicsModel * model;

  model = new akantu::StructuralMechanicsModel(beams);

  StructuralMaterial mat1;
  mat1.E=3e10;
  mat1.I=0.0025;
  mat1.A=0.01;


  model->addMaterial(mat1);

  StructuralMaterial mat2 ;
  mat2.E=3e10;
  mat2.I=0.003125;
  mat2.A=0.01;

  model->addMaterial(mat2);

  /* -------------------------------------------------------------------------- */
  // Defining the forces

  model->initModel();

  model->initVectors();

  UInt nb_element = beams.getNbElement(type);
  for (unsigned int i = 0; i < nb_element; ++i) {
    model->getElementMaterial(type)(i,0) = beams.getUIntData(type, "tag_0")(i,0) - 1;
  }


  Vector<Real> & forces = model->getForce();
  Vector<Real> & displacement = model->getDisplacement();
  Vector<bool> & boundary = model->getBoundary();
  const Vector<Real> & N_M  = model->getStress(_bernoulli_beam_2);

  Vector<UInt> & element_material = model->getElementMaterial(_bernoulli_beam_2);

  forces.clear();
  displacement.clear();


  model->computeForcesFromFunction(lin_load, akantu::_bft_forces);

  /* -------------------------------------------------------------------------- */
  // Defining the boundary conditions

  boundary(0,0) = true;
  boundary(0,1) = true;
  boundary(3,0) = true;
  boundary(3,1) = true;
  boundary(4,0) = true;
  boundary(4,1) = true;
  boundary(4,2) = true;
  boundary(5,0) = true;
  boundary(5,1) = true;
  boundary(5,2) = true;
  //boundary(2,1) = true;
  //boundary(2,0) = true;
  //boundary(1,1) = true;
  //boundary(1,0) = true;


  /* -------------------------------------------------------------------------- */
  // Solve

  model->initImplicitSolver();

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

  /* -------------------------------------------------------------------------- */
  // Post-Processing

  model->computeStressOnQuad();

  model->getStiffnessMatrix().saveMatrix("Ka.mtx");
  std::cout<< " x1 = " << displacement(1,2) << std::endl;
  std::cout<< " x2 = " << displacement(2,2) << std::endl;


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
template <> UInt paraviewType<_bernoulli_beam_2>(){ return BEAM2; };
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
  dumper.AddElemDataField(model.getStress(_bernoulli_beam_2).values,
			  2, "stress");
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
