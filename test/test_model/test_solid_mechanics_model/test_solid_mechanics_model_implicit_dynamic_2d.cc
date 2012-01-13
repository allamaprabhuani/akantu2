/**
 * @file   test_solid_mechanics_model_implicit_dynamic_2d.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Apr 29 11:32:25 2011
 *
 * @brief  test of the dynamic implicit code
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
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "static_communicator.hh"
#include "distributed_synchronizer.hh"
#include "mesh_partition_scotch.hh"

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

using namespace akantu;

/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"

void paraviewInit(iohelper::Dumper & dumper, const SolidMechanicsModel & model);
void paraviewDump(iohelper::Dumper & dumper);
#endif

/* -------------------------------------------------------------------------- */
const Real F = 0.5e4;
#define bar_length 10.
#define bar_height 1.
#define bar_depth  1.
const ElementType TYPE = _triangle_3;

UInt spatial_dimension = 2;
Real time_step = 1e-4;


Real analytical_solution(Real time) {
  return 1./pow(M_PI, 4) * ((1. - cos(M_PI*M_PI*time)) + (1. - cos(3*3*M_PI*M_PI*time))/81. + (1. - cos(5*5*M_PI*M_PI*time))/625.);
}

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  debug::setDebugLevel(dblWarning);
  initialize(argc, argv);

  Mesh mesh(spatial_dimension);

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm->getNbProc();
  Int prank = comm->whoAmI();

  MeshPartition * partition = NULL;
  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("beam_2d_lin.msh", mesh);
    
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->reorder();
    partition->partitionate(psize);
  }
  
  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);
  model->initParallel(partition);
  

  //  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();

  /// model initialization
  model->initVectors();

  model->initModel();
  model->readMaterials("material_implicit_dynamic.dat");
  model->initMaterials();

  model->initImplicit(true);

  // boundary conditions
  const Vector<Real> & position    = mesh.getNodes();
  Vector<bool> & boundary    = model->getBoundary();
  Vector<Real> & force       = model->getForce();
  Vector<Real> & displacment = model->getDisplacement();

  //initial conditions
  model->getForce().clear();
  model->getDisplacement().clear();
  model->getVelocity().clear();
  model->getAcceleration().clear();

  // MeshUtils::buildFacets(mesh);
  // MeshUtils::buildSurfaceID(mesh);

  // CSR<UInt> surface_nodes;
  // MeshUtils::buildNodesPerSurface(mesh, surface_nodes);

  UInt node_to_print = -1;
  bool print_node = false;

  // for (UInt s = 0; s < surface_nodes.getNbRows(); ++s) {
  //   CSR<UInt>::iterator snode = surface_nodes.begin(s);
  //   for(; snode != surface_nodes.end(s); ++snode) {
  //    UInt n = *snode;
  Vector<UInt> node_to_displace;
  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    Real x = position(n, 0);
    Real y = position(n, 1);
    Real z = 0;
    if(spatial_dimension == 3)
      z = position(n, 2);

    if(Math::are_float_equal(x, 0.) &&
       Math::are_float_equal(y, bar_height / 2.)) {
      boundary(n,0) = true;
      boundary(n,1) = true;
      if(spatial_dimension == 3 && Math::are_float_equal(z, bar_depth / 2.))
	boundary(n,2) = true;
    }

    if(Math::are_float_equal(x, bar_length) &&
       Math::are_float_equal(y, bar_height / 2.)) {
      boundary(n,1) = true;
      if(spatial_dimension == 3 && Math::are_float_equal(z, bar_depth / 2.))
	boundary(n,2) = true;
    }

    if(Math::are_float_equal(x, bar_length / 2.) &&
       Math::are_float_equal(y, bar_height / 2.)) {

      if(spatial_dimension < 3 || (spatial_dimension == 3 && Math::are_float_equal(z, bar_depth / 2.))) {
	force(n,1) = F;
	if(mesh.isLocalOrMasterNode(n)) {
	  print_node = true;
	  node_to_print = n;
	  std::cout << "I, proc " << prank +1 << " handle the print of node " << n 
		    << "(" << x << ", "<< y << ", " << z << ")" << std::endl;
	}
      }
    }
  }
  //  }

  model->setTimeStep(time_step);

  model->updateResidual();

  std::stringstream out;
  out << "position-" << TYPE << "_" << std::scientific << time_step << ".csv";

#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, *model);
#endif

  std::ofstream pos;
  if(print_node) {
    pos.open(out.str().c_str());
    if(!pos.good()) {
      std::cerr << "Cannot open file " << out.str() <<std::endl;
      exit(EXIT_FAILURE);
    }
    pos << "id,time,position,solution" << std::endl;
  }

  Real time = 0;
  UInt count = 0;
  //  UInt print_freq = 1;
  Real error;

  model->assembleStiffnessMatrix();
  model->assembleMass();

  // model->assembleMassLumped();
  // Vector<Real> lumped_mass(0,spatial_dimension);
  // model->getMassMatrix().lump(lumped_mass);

  // debug::setDebugLevel(dblTest);
  // std::cout << model->getMass() << lumped_mass;
  // debug::setDebugLevel(dblWarning);

  model->getMassMatrix().saveMatrix("M.mtx");
  model->getStiffnessMatrix().saveMatrix("K.mtx");

  /// time loop
  for (UInt s = 1; time < 0.62; ++s) {
    model->implicitPred();
    /// convergence loop
    do {
      if(count > 0 && prank == 0)
    	std::cout << "passing step " << s << " " << s * time_step << "s - " << std::setw(4) << count << " : " << std::scientific << error << "\r" << std::flush;
      model->updateResidual();
      model->solveDynamic();
      model->implicitCorr();
      count++;
    } while(!model->testConvergenceIncrement(1e-12, error) && (count < 1000));
    if(prank == 0) std::cout << "passing step " << s << " " << s * time_step << "s - " << std::setw(4) << count << " : " << std::scientific << error << std::endl;
    count = 0;

    // if(s % print_freq == 0) {
    //   std::cout << "passing step " << s << "/" << max_steps << " " << s * time_step << "s - " << count / print_freq << std::endl;
    //   count = 0;
    // }

    if(print_node) pos << s << "," << s * time_step << "," << displacment(node_to_print,  1) << "," << analytical_solution(s*time_step) << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) paraviewDump(dumper);
#endif

    time += time_step;
  }

  if(print_node) pos.close();

  delete model;
  finalize();

  return EXIT_SUCCESS;
}



/* -------------------------------------------------------------------------- */
/* iohelper::Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
template <ElementType type> iohelper::ElemType paraviewType();
template <> iohelper::ElemType paraviewType<_segment_2>()      { return iohelper::LINE1;     };
template <> iohelper::ElemType paraviewType<_segment_3>()      { return iohelper::LINE2;     };
template <> iohelper::ElemType paraviewType<_triangle_3>()     { return iohelper::TRIANGLE1; };
template <> iohelper::ElemType paraviewType<_triangle_6>()     { return iohelper::TRIANGLE2; };
template <> iohelper::ElemType paraviewType<_quadrangle_4>()   { return iohelper::QUAD1;     };
template <> iohelper::ElemType paraviewType<_tetrahedron_4>()  { return iohelper::TETRA1;    };
template <> iohelper::ElemType paraviewType<_tetrahedron_10>() { return iohelper::TETRA2;    };
template <> iohelper::ElemType paraviewType<_hexahedron_8>()   { return iohelper::HEX1;      };

/* -------------------------------------------------------------------------- */
void paraviewInit(iohelper::Dumper & dumper, const SolidMechanicsModel & model) {
  UInt spatial_dimension = ElementClass<TYPE>::getSpatialDimension();
  UInt nb_nodes   = model.getFEM().getMesh().getNbNodes();
  UInt nb_element = model.getFEM().getMesh().getNbElement(TYPE);

  std::stringstream filename; filename << "dynamic_implicit_beam_" << TYPE;

  dumper.SetMode(iohelper::TEXT);
  dumper.SetParallelContext(StaticCommunicator::getStaticCommunicator()->whoAmI(), StaticCommunicator::getStaticCommunicator()->getNbProc());
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, filename.str().c_str());
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(TYPE).values,
			 paraviewType<TYPE>(), nb_element, iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getAcceleration().values,
			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model.getResidual().values,
			  spatial_dimension, "residual");
  dumper.AddNodeDataField(model.getForce().values,
			  spatial_dimension, "applied_force");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(TYPE).values,
   			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(TYPE).values,
   			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(iohelper::Dumper & dumper) {
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */

#endif
