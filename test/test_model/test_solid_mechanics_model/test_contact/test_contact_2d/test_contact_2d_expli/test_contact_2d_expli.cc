/**
 * @file   test_contact_2d_expli.cc
 *
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 *
 * @date   Mon Jan 17 12:17:38 2011
 *
 * @brief  test explicit DCR contact algorithm for 2d
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
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "contact.hh"
#include "contact_2d_explicit.hh"
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
#include "io_helper.hh"

#endif //AKANTU_USE_IOHELPER

#define NORMAL_PRESSURE -1.e6

using namespace akantu;

//static void reduceGap(const SolidMechanicsModel & model, const Real threshold, const Real gap);
static void setBoundaryConditions(SolidMechanicsModel & model);
static void reduceVelocities(const SolidMechanicsModel & model, const Real ratio);
static void initParaview(SolidMechanicsModel & model);

#ifdef AKANTU_USE_IOHELPER
iohelper::DumperParaview dumper;
#endif

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);

  UInt spatial_dimension = 2;
  UInt max_steps = 30000;
  Real time_factor = 0.2;

  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("squares.msh", mesh);

  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);

  /// get two squares closer
  // reduceGap(*model, 0.05, 1.e-6);

  UInt nb_nodes = model->getFEEngine().getMesh().getNbNodes();
  UInt nb_elements = model->getFEEngine().getMesh().getNbElement(_triangle_3);

  /// model initialization
  model->initModel();
  model->initArrays();

  model->readMaterials("materials.dat");
  model->initMaterials();
  model->initExplicit();
  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMassLumped();

  /// set vectors to zero
  memset(model->getForce().storage(),        0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getVelocity().storage(),     0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getAcceleration().storage(), 0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getDisplacement().storage(), 0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getResidual().storage(), 0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getMaterial(0).getStrain(_triangle_3).storage(), 0,
	 spatial_dimension*spatial_dimension*nb_elements*sizeof(Real));
  memset(model->getMaterial(0).getStress(_triangle_3).storage(), 0,
	 spatial_dimension*spatial_dimension*nb_elements*sizeof(Real));

  /// Paraview Helper
  #ifdef AKANTU_USE_IOHELPER
  initParaview(*model);
  #endif //AKANTU_USE_IOHELPER

  Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

  /// set boundary conditions
  setBoundaryConditions(*model);

  /// define and initialize contact
  Contact * contact = Contact::newContact(*model,
					  _ct_2d_expli,
					  _cst_2d_expli,
					  _cnst_2d_grid);

  Contact2dExplicit * my_contact = dynamic_cast<Contact2dExplicit *>(contact);

  my_contact->initContact(true);
  my_contact->setFrictionCoefficient(0.);
  my_contact->initNeighborStructure();
  my_contact->initSearch();

  for (UInt s = 0; s < max_steps; ++s) {

    model->explicitPred();

    model->updateCurrentPosition();

    my_contact->solveContact();

    model->updateResidual();

    model->updateAcceleration();
    model->explicitCorr();

#ifdef AKANTU_USE_IOHELPER
    if(s % 200 == 0)
      dumper.Dump();
#endif

    if(s%100 == 0 && s>499)
      reduceVelocities(*model, 0.95);

    if(s % 500 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  delete my_contact;
  delete model;
  finalize();
  return EXIT_SUCCESS;
}



// /* -------------------------------------------------------------------------- */
// static void reduceGap(const SolidMechanicsModel & model, const Real threshold, const Real gap) {

//   UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
//   Real * coord = model.getFEEngine().getMesh().getNodes().storage();
//   Real y_top = HUGE_VAL, y_bot = -HUGE_VAL;

//   for (UInt n = 0; n < nb_nodes; ++n) {
//     if (coord[2*n+1] > threshold) {
//       if(coord[2*n+1] < y_top)
// 	y_top = coord[2*n+1];
//     }
//     else {
//       if (coord[2*n+1] > y_bot)
// 	y_bot = coord[2*n+1];
//     }
//   }

//   Real delta = y_top - y_bot - gap;
//   /// move all nodes belonging to the top cube
//   for (UInt n = 0; n < nb_nodes; ++n) {
//     if (coord[2*n+1] > threshold)
//       coord[2*n+1] -= delta;
//   }
// }

class MyStressFunctor : public SolidMechanicsModel::SurfaceLoadFunctor {
public:
  MyStressFunctor(Real y_max) : y_max(y_max) {};

  inline void stress(const Vector<Real> & position,
		     Matrix & stress,
		     __attribute__ ((unused)) const Vector<Real> & normal,
		     __attribute__ ((unused)) Surface surface_id) {
    stress.clear();
    if(position(1) > y_max - 1.e-5)
      stress(1,1) = NORMAL_PRESSURE;
  }
private:
  Real y_max;
};

/* -------------------------------------------------------------------------- */
static void setBoundaryConditions(SolidMechanicsModel & model) {

  UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  Real * coord = model.getFEEngine().getMesh().getNodes().storage();
  Real y_min = std::numeric_limits<Real>::max(), y_max = -std::numeric_limits<Real>::max();
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (coord[2*n+1] > y_max)
      y_max = coord[2*n+1];
    if (coord[2*n+1] < y_min)
      y_min = coord[2*n+1];
  }

  FEEngine & b_fem = model.getFEEngineBoundary();
  b_fem.initShapeFunctions();
  b_fem.computeNormalsOnControlPoints();
  bool * id = model.getBlockedDOFs().storage();
  memset(id, 0, 2*nb_nodes*sizeof(bool));
  std::cout << "Nodes ";
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (coord[2*i+1] < y_min + 1.e-5) {
      id[2*i+1] = true;
      std::cout << " " << i << " ";
    }
  }
  std::cout << "are blocked" << std::endl;

  MyStressFunctor func(y_max);
  model.computeForcesFromFunction(func, _bft_stress);
}

/* -------------------------------------------------------------------------- */
/// artificial damping of velocities in order to reach a global static equilibrium
static void reduceVelocities(const SolidMechanicsModel & model, const Real ratio)
{
  UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  Real * velocities = model.getVelocity().storage();

  if(ratio>1.) {
    fprintf(stderr,"**error** in Reduce_Velocities ratio bigger than 1!\n");
    exit(-1);
  }

  for(UInt i =0; i<nb_nodes; i++) {
    velocities[2*i] *= ratio;
    velocities[2*i+1] *= ratio;
  }
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
static void initParaview(SolidMechanicsModel & model)
{
  UInt spatial_dimension = model.getSpatialDimension();
  UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  UInt nb_elements = model.getFEEngine().getMesh().getNbElement(_triangle_3);

  dumper.SetMode(iohelper::TEXT);
  dumper.SetPoints(model.getFEEngine().getMesh().getNodes().storage(),
		   spatial_dimension, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model.getFEEngine().getMesh().getConnectivity(_triangle_3).storage(),
			 iohelper::TRIANGLE1, nb_elements, iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().storage(),
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getVelocity().storage(),
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getResidual().storage(),
			  spatial_dimension, "force");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(_triangle_3).storage(),
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStress(_triangle_3).storage(),
			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}
#endif
