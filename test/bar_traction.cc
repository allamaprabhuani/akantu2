/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

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
  akantu::UInt max_steps = 50000;
  akantu::Real time_factor = 0.2;

  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("bar.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  akantu::UInt nb_nodes = model->getFEM().getNbNodes();
  akantu::UInt nb_element = model->getFEM().getNbElement(akantu::_triangle_1);

  /// model initialization
  model->initVectors();

  /// set vectors to 0
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));

#ifdef AKANTU_USE_IOHELPER
  /// set to 0 only for the first paraview dump
  memset(model->getResidual().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getStrain(akantu::_triangle_1).values, 0,
	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
  memset(model->getStress(akantu::_triangle_1).values, 0,
	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
#endif //AKANTU_USE_IOHELPER

  model->readMaterials("");
  model->initMaterials();
  model->initModel();

  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMass();

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


#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(akantu::_triangle_1).values,
			 TRIANGLE1, nb_element, C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddElemDataField(model->getStrain(akantu::_triangle_1).values,
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model->getStress(akantu::_triangle_1).values,
			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
#endif //AKANTU_USE_IOHELPER

  dumper.Dump();

  model->setPotentialEnergyFlagOn();
  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    // epot = model->getPotentialEnergy();
    // ekin = model->getKineticEnergy();
    // std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
    // 	      << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  return EXIT_SUCCESS;
}
