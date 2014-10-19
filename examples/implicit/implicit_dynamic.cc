/**
 * @file   implicit_dynamic.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Feb 24 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  This code refers to the implicit dynamic example from the user manual
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model.hh"

#include <iostream>
using namespace akantu;

/* -------------------------------------------------------------------------- */
const Real bar_length = 10.;
const Real bar_height = 1.;
const Real bar_depth  = 1.;
const Real F   = 5e3;
const Real L   = bar_length;
const Real I   = bar_depth * bar_height * bar_height * bar_height / 12.;
const Real E   = 12e7;
const Real rho = 1000;
const Real m   = rho * bar_height * bar_depth;

static Real w(UInt n) {
  return n*n*M_PI*M_PI/(L*L)*sqrt(E*I/m);
}

static Real analytical_solution(Real time) {
  return 2*F*L*L*L/(pow(M_PI, 4)*E*I) * ((1. - cos(w(1)*time)) + (1. - cos(w(3)*time))/81. + (1. - cos(w(5)*time))/625.);
}

const UInt spatial_dimension = 2;
const Real time_step = 1e-4;
const Real max_time  = 0.62;
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  initialize("material_dynamic.dat", argc, argv);

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  MeshPartition * partition = NULL;
  if(prank == 0) {
    mesh.read("beam.msh");
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  SolidMechanicsModel model(mesh);
  model.initParallel(partition);
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_implicit_dynamic));
  Material &mat = model.getMaterial(0);
  mat.setParam("E",   E);
  mat.setParam("rho", rho);

  model.getMassMatrix().saveMatrix("M.mtx");

  Array<Real> & force          = model.getForce();
  Array<Real> & displacment    = model.getDisplacement();

  // boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "blocked");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "blocked");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "roller");

  const Array<UInt> & trac_nodes = mesh.getElementGroup("traction").getNodes();

  bool dump_node = false;
  if(trac_nodes.getSize() > 0 && mesh.isLocalOrMasterNode(trac_nodes(0))) {
    force(trac_nodes(0), 1) = F;
    dump_node = true;
  }

  // output setup
  std::ofstream pos;
  pos.open("position.csv");
  if(!pos.good()) AKANTU_DEBUG_ERROR("Cannot open file \"position.csv\"");

  pos << "id,time,position,solution" << std::endl;

  model.setBaseName("dynamic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.dump();


  model.setTimeStep(time_step);

  /// time loop
  Real time = 0.;
  for (UInt s = 1; time < max_time; ++s, time += time_step) {
    if(prank == 0)
      std::cout << s << "\r" << std::flush;

    model.solveStep<_scm_newton_raphson_tangent_modified, _scc_increment>(1e-12, 100);

    if(dump_node)
      pos << s << ","
          << time << ","
          << displacment(trac_nodes(0),  1) << ","
          << analytical_solution(s*time_step) << std::endl;

    if(s % 100 == 0)
      model.dump();
  }

  std::cout << std::endl;
  pos.close();

  finalize();

  return EXIT_SUCCESS;
}
