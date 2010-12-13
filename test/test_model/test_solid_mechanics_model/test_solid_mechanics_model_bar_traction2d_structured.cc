/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * \<insert license here\>
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
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

int main(int argc, char *argv[])
{
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 10000;
  akantu::Real time_factor = 0.2;

  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;

  mesh_io.read("bar_structured1.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  /// model initialization
  model->initVectors();

  /// set vectors to 0
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));


  model->readMaterials("material.dat");
  model->initMaterials();

  model->initModel();

  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMassLumped();


  /// boundary conditions
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] >= 9)
      model->getDisplacement().values[spatial_dimension*i] = (model->getFEM().getMesh().getNodes().values[spatial_dimension*i] - 9) / 100. ;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] <= eps)
	model->getBoundary().values[spatial_dimension*i] = true;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] <= eps ||
       model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] >= 1 - eps ) {
      model->getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }

  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;


  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    epot = model->getPotentialEnergy();
    ekin = model->getKineticEnergy();

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;
  }

  energy.close();

  return EXIT_SUCCESS;
}
