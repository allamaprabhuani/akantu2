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
#include "communicator.hh"
#include "mesh_partition_scotch.hh"

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

using namespace akantu;

/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"

void paraviewInit(Dumper & dumper, const SolidMechanicsModel & model);
void paraviewDump(Dumper & dumper);
#endif

/* -------------------------------------------------------------------------- */

#define bar_length 10
#define bar_height 1
#define TYPE _tetrahedron_10

UInt spatial_dimension = 3;
Real time_step = 1e-4;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  debug::setDebugLevel(dblWarning);
  initialize(&argc, &argv);


  UInt max_steps = 10000000;


  Mesh mesh(spatial_dimension);

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm->getNbProc();
  Int prank = comm->whoAmI();

  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  Communicator * communicator;
  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("beam_3d_quad.msh", mesh);

    MeshPartition * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->reorder();
    partition->partitionate(psize);
    communicator = Communicator::createCommunicatorDistributeMesh(mesh, partition);
    delete partition;
  } else {
    communicator = Communicator::createCommunicatorDistributeMesh(mesh, NULL);
  }

  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);

  //  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();

  /// model initialization
  model->initVectors();

  model->initModel();
  model->readMaterials("material.dat");
  model->initMaterials();
  model->registerSynchronizer(*communicator);

  model->initImplicit(true);
  model->assembleMass();

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

  MeshUtils::buildFacets(mesh);
  MeshUtils::buildSurfaceID(mesh);

  CSR<UInt> surface_nodes;
  MeshUtils::buildNodesPerSurface(mesh, surface_nodes);

  for (UInt s = 0; s < surface_nodes.getNbRows(); ++s) {
    CSR<UInt>::iterator snode = surface_nodes.begin(s);
    for(; snode != surface_nodes.end(s); ++snode) {
     UInt n = *snode;
     Real x = position(n, 0);
     Real y = position(n, 1);
     Real z = -1;
     if(spatial_dimension == 3)
       z = position(n, 2);

     if(Math::are_float_equal(x, 0.) &&
	Math::are_float_equal(y, 0.)) {
       boundary(n,0) = true;
       boundary(n,1) = true;
       if(spatial_dimension == 3)
	 boundary(n,2) = true;
     }

     if(Math::are_float_equal(x, bar_length) &&
	Math::are_float_equal(y, 0.)) {
       boundary(n,1) = true;
       if(spatial_dimension == 3)
	 boundary(n,2) = true;
     }

     if(Math::are_float_equal(x, bar_length / 2.) &&
	Math::are_float_equal(y, bar_height)) {
       force(n,1) = 10000;
     }
    }
  }

  model->setTimeStep(time_step);

  model->updateResidual();
  //  model->initialAcceleration();

  std::stringstream out;
  out << "position_" << std::scientific << time_step << ".csv";

  DumperParaview dumper;
  paraviewInit(dumper, *model);

  std::ofstream pos;
  pos.open(out.str().c_str());
  pos << "id,time,position" << std::endl;

  Real time = 0;
  UInt count = 0;
  //  UInt print_freq = 1;
  Real error;

  model->getMassMatrix().saveMatrix("M.mtx");
  model->assembleStiffnessMatrix();

  /// time loop
  for (UInt s = 1; s <= max_steps && time < 0.05; ++s) {
    model->implicitPred();

    /// convergence loop
    do {
      if(count > 0) 
	std::cout << "passing step " << s << "/" << max_steps << " " << s * time_step << "s - " << std::setw(4) << count << " : " << std::scientific << error << "\r" << std::flush;
      model->updateResidual();
      model->solveDynamic();
      model->implicitCorr();
      count++;

    } while(!model->testConvergenceIncrement(1e-19, error) && (count < 10000));
    std::cout << "passing step " << s << "/" << max_steps << " " << s * time_step << "s - " << std::setw(4) << count << " : " << std::scientific << error << std::endl;
    count = 0;

    // if(s % print_freq == 0) {
    //   std::cout << "passing step " << s << "/" << max_steps << " " << s * time_step << "s - " << count / print_freq << std::endl;
    //   count = 0;
    // }

    pos << s << "," << s * time_step << "," << displacment(3,1) << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) paraviewDump(dumper);
#endif

    time += time_step;
  }

  delete model;
  finalize();

  return EXIT_SUCCESS;
}



/* -------------------------------------------------------------------------- */
/* Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
template <ElementType type> UInt paraviewType();
template <> UInt paraviewType<_segment_2>()      { return LINE1;     };
template <> UInt paraviewType<_segment_3>()      { return LINE2;     };
template <> UInt paraviewType<_triangle_3>()     { return TRIANGLE1; };
template <> UInt paraviewType<_triangle_6>()     { return TRIANGLE2; };
template <> UInt paraviewType<_quadrangle_4>()   { return QUAD1;     };
template <> UInt paraviewType<_tetrahedron_4>()  { return TETRA1;    };
template <> UInt paraviewType<_tetrahedron_10>() { return TETRA2;    };
template <> UInt paraviewType<_hexahedron_8>()   { return HEX1;      };

/* -------------------------------------------------------------------------- */
void paraviewInit(Dumper & dumper, const SolidMechanicsModel & model) {
  UInt spatial_dimension = ElementClass<TYPE>::getSpatialDimension();
  UInt nb_nodes   = model.getFEM().getMesh().getNbNodes();
  UInt nb_element = model.getFEM().getMesh().getNbElement(TYPE);

  std::stringstream filename; filename << "dynamic_implicit_beam_" << TYPE;

  dumper.SetMode(TEXT);
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, filename.str().c_str());
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(TYPE).values,
			 paraviewType<TYPE>(), nb_element, C_MODE);
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
void paraviewDump(Dumper & dumper) {
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */

#endif
