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
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "fem.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

static void trac(__attribute__ ((unused)) double * position,double * stress){
  memset(stress,0,sizeof(Real)*4);
  stress[0] = 1000;
  stress[3] = 1000;
}

int main(int argc, char *argv[])
{
  UInt max_steps = 10000;
  Real epot, ekin;

  Mesh mesh(2);
  MeshIOMSH mesh_io;
  mesh_io.read("circle2.msh", mesh);

  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);

  /// model initialization
  model->initVectors();
  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0, 2*nb_nodes*sizeof(Real));
  memset(model->getVelocity().values,     0, 2*nb_nodes*sizeof(Real));
  memset(model->getAcceleration().values, 0, 2*nb_nodes*sizeof(Real));
  memset(model->getDisplacement().values, 0, 2*nb_nodes*sizeof(Real));
  memset(model->getResidual().values,     0, 2*nb_nodes*sizeof(Real));
  memset(model->getMass().values,     1, nb_nodes*sizeof(Real));

  model->readCustomMaterial<MaterialDamage>("material.dat","DAMAGE");

  model->initMaterials();
  model->initModel();

  Real time_step = model->getStableTimeStep();
  model->setTimeStep(time_step/10.);

  model->assembleMassLumped();

  std::cout << *model << std::endl;

  FEM & fem_boundary = model->getFEMBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnQuadPoints();
  model->computeForcesFromFunction(trac, akantu::_bft_forces);

#ifdef AKANTU_USE_IOHELPER
  model->updateResidual();

  DumperParaview dumper;
  dumper.SetMode(BASE64);

  dumper.SetPoints(model->getFEM().getMesh().getNodes().values, 2, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(_triangle_6).values,
			 TRIANGLE2, model->getFEM().getMesh().getNbElement(_triangle_6), C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values, 2, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values, 2, "velocity");
  dumper.AddNodeDataField(model->getForce().values, 2, "force");
  dumper.AddNodeDataField(model->getMass().values, 1, "Mass");
  dumper.AddNodeDataField(model->getResidual().values, 2, "residual");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(_triangle_6).values, 4, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(_triangle_6).values, 4, "stress");
  MaterialDamage & mat = dynamic_cast<MaterialDamage&>(model->getMaterial(0));
  Real * dam = mat.getDamage(_triangle_6).values;
  dumper.AddElemDataField(dam, 1, "damage");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("force", 1);
  dumper.SetEmbeddedValue("residual", 1);
  dumper.SetEmbeddedValue("velocity", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  for(UInt s = 0; s < max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    epot = model->getPotentialEnergy();
    ekin = model->getKineticEnergy();

    std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  return EXIT_SUCCESS;
}



