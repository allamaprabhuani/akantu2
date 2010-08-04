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
  akantu::UInt max_steps = 10000;
  akantu::Real epot, ekin, epot_first;
  bool first = true;

  akantu::Mesh mesh(1);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("line.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(1, mesh);

  /// model initialization
  model->initVectors();
  model->readMaterials("");
  model->initMaterials();
  model->initModel();

  akantu::Real time_step = model->getStableTimeStep();
  model->setTimeStep(time_step);

  model->assembleMass();

  akantu::UInt nb_nodes = model->getFEM().getNbNodes();

  /// boundary conditions
  memset(model->getForce().values,        0, nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0, nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0, nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0, nb_nodes*sizeof(akantu::Real));
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    model->getDisplacement().values[i] = model->getFEM().getMesh().getNodes().values[i] / 100.;
  }
  //  model->getDisplacement().values[1] = 0.1;
  model->getBoundary().values[0] = true;


#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);

  dumper.SetPoints(model->getFEM().getMesh().getNodes().values, 1, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(akantu::_line_1).values,
			 LINE1, model->getFEM().getNbElement(akantu::_line_1), C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values, 1, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values, 1, "velocity");
  dumper.AddNodeDataField(model->getResidual().values, 1, "force");
  dumper.AddElemDataField(model->getStrain(akantu::_line_1).values, 1, "strain");
  dumper.AddElemDataField(model->getStress(akantu::_line_1).values, 1, "stress");
  //  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
#endif //AKANTU_USE_IOHELPER

  model->setPotentialEnergyFlagOn();
  for(akantu::UInt s = 0; s < max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();


    epot = model->getPotentialEnergy();
    ekin = model->getKineticEnergy();
    if(first) {
      epot_first = epot;
      first = false;
    }

    std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
	      << " " << model->getDisplacement().values[1]
	      << " " << model->getVelocity().values[1]
	      << " " << model->getResidual().values[1]
	      << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  return EXIT_SUCCESS;
}
