/**
 * @file   test_solid_mechanics_model_cube3d.cc
 * @author Guillaume ANCIAUX <anciaux@lsmscluster1.epfl.ch>
 * @date   Tue Aug 17 11:31:22 2010
 *
 * @brief  test of the class SolidMechanicsModel on the 3d cube
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

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
  akantu::Real epot, ekin;

  akantu::Mesh mesh(3);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("cube.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  /// model initialization
  model->initVectors();
  /// initialize the vectors
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0, 3*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0, 3*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0, 3*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0, 3*nb_nodes*sizeof(akantu::Real));

  model->readMaterials("material.dat");
  model->initMaterials();
  model->initModel();

  akantu::Real time_step = model->getStableTimeStep();
  model->setTimeStep(time_step/10.);

  model->assembleMassLumped();

  std::cout << *model << std::endl;


  /// boundary conditions
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    model->getDisplacement().values[3*i] = model->getFEM().getMesh().getNodes().values[3*i] / 100.;

    if(model->getFEM().getMesh().getNodes().values[3*i] <= eps) {
      model->getBoundary().values[3*i    ] = true;
      if(model->getFEM().getMesh().getNodes().values[3*i + 1] <= eps)
	model->getBoundary().values[3*i + 1] = true;
    }
    if(model->getFEM().getMesh().getNodes().values[3*i + 1] <= eps) {
      model->getBoundary().values[3*i + 1] = true;
    }

  }
  //  model->getDisplacement().values[1] = 0.1;


#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  //  dumper.SetMode(TEXT);

  dumper.SetPoints(model->getFEM().getMesh().getNodes().values, 3, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(akantu::_tetrahedron_6).values,
			 TETRA1, model->getFEM().getMesh().getNbElement(akantu::_tetrahedron_6), C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values, 3, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values, 3, "velocity");
  dumper.AddNodeDataField(model->getMass().values, 1, "mass");
  dumper.AddNodeDataField(model->getResidual().values, 3, "force");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(akantu::_tetrahedron_6).values, 9, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(akantu::_tetrahedron_6).values, 9, "stress");
  dumper.SetPrefix("paraview/");
  dumper.Init();
#endif //AKANTU_USE_IOHELPER

  for(akantu::UInt s = 0; s < max_steps; ++s) {
    model->explicitPred();
    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();


    epot = model->getPotentialEnergy();
    ekin = model->getKineticEnergy();

    std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  return EXIT_SUCCESS;
}
