/**
 * @file   test_solid_mechanics_model_cube3d.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Tue Aug 17 11:31:22 2010
 *
 * @brief  test of the class SolidMechanicsModel on the 3d cube
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
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER


int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::UInt max_steps = 10000;
  akantu::Real epot, ekin;

#ifdef AKANTU_USE_IOHELPER
  akantu::ElementType type = akantu::_tetrahedron_4;
  iohelper::ElemType paratype = iohelper::TETRA1;
#endif //AKANTU_USE_IOHELPER

  akantu::Mesh mesh(3);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("cube1.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  /// model initialization
  model->initVectors();
  /// initialize the vectors
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0, 3*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0, 3*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0, 3*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0, 3*nb_nodes*sizeof(akantu::Real));

  model->initExplicit();
  model->initModel();
  model->readMaterials("material.dat");
  model->initMaterials();

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
    }

    if(model->getFEM().getMesh().getNodes().values[3*i + 1] <= eps) {
      model->getBoundary().values[3*i + 1] = true;
    }

  }
  //  model->getDisplacement().values[1] = 0.1;


#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  //  dumper.SetMode(iohelper::TEXT);

  dumper.SetPoints(model->getFEM().getMesh().getNodes().values, 3, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 paratype, model->getFEM().getMesh().getNbElement(type), iohelper::C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values, 3, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values, 3, "velocity");
  dumper.AddNodeDataField(model->getMass().values, 1, "mass");
  dumper.AddNodeDataField(model->getResidual().values, 3, "force");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(type).values, 9, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(type).values, 9, "stress");
  dumper.SetPrefix("paraview/");
  dumper.Init();
#endif //AKANTU_USE_IOHELPER

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 0; s < max_steps; ++s) {
    model->explicitPred();
    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();


    epot = model->getPotentialEnergy();
    ekin = model->getKineticEnergy();

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  energy.close();

  akantu::finalize();
  return EXIT_SUCCESS;
}
