/**
 * @file   damage_non_local.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Sep 30 10:06:56 2011
 *
 * @brief
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

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "static_communicator.hh"
#include "mesh_partition_scotch.hh"
#include "io_helper_tools.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

const UInt spatial_dimension = 2;
const Real bar_length = 10.;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  debug::setDebugLevel(akantu::dblWarning);
  akantu::initialize(&argc, &argv);

  ElementType type = _triangle_6;
  UInt max_steps = 15000;

  Mesh mesh(spatial_dimension);

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm->getNbProc();
  Int prank = comm->whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("bar.msh", mesh);
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  akantu::SolidMechanicsModel model(mesh);
  model.initParallel(partition);

  model.initModel();
  model.initVectors();
  model.initExplicit();

  model.readMaterials("material.dat");
  model.initMaterials();

  /// set vectors to 0
  model.getForce().clear();
  model.getVelocity().clear();
  model.getAcceleration().clear();
  model.getDisplacement().clear();

  /// boundary conditions
  Real eps = 1e-16;
  const akantu::Vector<akantu::Real> & pos = mesh.getNodes();
  akantu::Vector<akantu::Real> & disp = model.getDisplacement();
  akantu::Vector<bool> & boun = model.getBoundary();

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    //    disp(i, 0) = pos(i,0) * 0.1;
    if(std::abs(pos(i, 0) - bar_length) <= eps) boun(i, 0) = true;
    if(std::abs(pos(i, 0) - 0) <= eps) boun(i, 0) = true;
  }

  model.assembleMassLumped();

  Real time_step = model.getStableTimeStep() * .8;
  if(prank == 0)
    std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  model.updateResidual();

  DumperParaview dumper;
  paraviewInit(dumper, model, type, "damage_nl_fine");
  //  paraviewInit(dumper, model, type, "damage");

  std::ofstream energy;
  energy.open("energy_damage_nl_fine.csv");
  //  energy.open("energy_damage.csv");

  energy << "id,rtime,epot,ekin,tot" << std::endl;

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    Real apply_disp = std::min(10e-4, (s-1)*time_step * 1e-6*1e6);
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if(std::abs(pos(i, 0) - bar_length) <= eps) disp(i, 0) = apply_disp;
      if(std::abs(pos(i, 0) - 0.) <= eps) disp(i, 0) = -apply_disp;
    }
    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    Real epot = model.getPotentialEnergy();
    Real ekin = model.getKineticEnergy();

    energy << s << "," << (s-1)*time_step << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

    if(s % 100 == 0) {
      std::cout << "." << std::flush;
      paraviewDump(dumper);
    }
  }
  std::cout << std::endl;
  energy.close();

  akantu::finalize();

  return EXIT_SUCCESS;
}
