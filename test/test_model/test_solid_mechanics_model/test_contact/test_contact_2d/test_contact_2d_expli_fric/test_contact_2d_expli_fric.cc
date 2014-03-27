/**
 * @file   test_contact_2d_expli_fric.cc
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
#include "io_helper.hh"
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
#include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

#define NORMAL_PRESSURE -1.e6
#define IMP_VEL 1.

using namespace akantu;

class Boundary {
public:
  Boundary(Mesh & mesh, SolidMechanicsModel & model);
  virtual ~Boundary();

  reduceGap(const Real threshold, const Real gap);
  setBoundaryConditions();

public:


public:
private:
  ///
  SolidMechanicsModel & model;

  ///
  Mesh & mesh;

  Real top_bounds[];
  Real
};



static void
static void setBoundaryConditions(SolidMechanicsModel & model);
void my_force(double * coord, double *T);
static void reduceVelocities(const SolidMechanicsModel & model, const Real ratio);
static void initParaview(SolidMechanicsModel & model);

Real y_min, y_max;
iohelper::DumperParaview dumper;

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

  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  UInt nb_elements = model->getFEM().getMesh().getNbElement(_triangle_3);

  /// model initialization
  model->initArrays();

  model->readMaterials("materials.dat");
  model->initMaterials();

  model->initModel();
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
  Boundary * my_boudary;
  setBoundaryConditions(*model);

  /// define and initialize contact
  Contact * my_contact = Contact::newContact(*model,
					     _ct_2d_expli,
					     _cst_2d_expli,
					     _cnst_2d_grid);

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

    if(s % 200 == 0)
      dumper.Dump();

    if(s%100 == 0 && s>499)
      reduceVelocities(*model, 0.95);

    if(s % 500 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  delete my_boudary;
  delete my_contact;
  delete model;
  finalize();
  return EXIT_SUCCESS;
}



/* -------------------------------------------------------------------------- */
static void reduceGap(const SolidMechanicsModel & model, const Real threshold, const Real gap) {

  UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  Real * coord = model.getFEM().getMesh().getNodes().storage();
  Real y_top = HUGE_VAL, y_bot = -HUGE_VAL;

  for (UInt n = 0; n < nb_nodes; ++n) {
    if (coord[2*n+1] > threshold) {
      if(coord[2*n+1] < y_top)
	y_top = coord[2*n+1];
    }
    else {
      if (coord[2*n+1] > y_bot)
	y_bot = coord[2*n+1];
    }
  }

  Real delta = y_top - y_bot - gap;
  /// move all nodes belonging to the top cube
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (coord[2*n+1] > threshold)
      coord[2*n+1] -= delta;
  }
}

/* -------------------------------------------------------------------------- */
static void setBoundaryConditions(SolidMechanicsModel & model) {

  UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  Real * coord = model.getFEM().getMesh().getNodes().storage();
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (coord[2*n+1] > y_max)
      y_max = coord[2*n+1];
    if (coord[2*n+1] < y_min)
      y_min = coord[2*n+1];
  }

  FEM & b_fem = model.getFEMBoundary();
  b_fem.initShapeFunctions();
  b_fem.computeNormalsOnQuadPoints();
  bool * id = model.getBoundary().storage();
  memset(id, 0, 2*nb_nodes*sizeof(bool));
  std::cout << "Nodes ";
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (coord[2*i+1] < y_min + 1.e-5) {
      id[2*i+1] = true;
      std::cout << " " << i << " ";
    }
  }
  std::cout << "are blocked" << std::endl;

  model.computeForcesFromFunction(my_force, _bft_stress);
}

/* -------------------------------------------------------------------------- */
void my_force(double * coord, double *T) {

  memset(T, 0, 4*sizeof(double));
  if(*(coord+1) > y_max-1.e-5)
    T[3] = NORMAL_PRESSURE;
}

/* -------------------------------------------------------------------------- */
/// artificial damping of velocities in order to reach a global static equilibrium
static void reduceVelocities(const SolidMechanicsModel & model, const Real ratio)
{
  UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
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
static void initParaview(SolidMechanicsModel & model)
{
  UInt spatial_dimension = model.getSpatialDimension();
  UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  UInt nb_elements = model.getFEM().getMesh().getNbElement(_triangle_3);

  dumper.SetMode(iohelper::TEXT);
  dumper.SetPoints(model.getFEM().getMesh().getNodes().storage(),
		   spatial_dimension, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(_triangle_3).storage(),
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
