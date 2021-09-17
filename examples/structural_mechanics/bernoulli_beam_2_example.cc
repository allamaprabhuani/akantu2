/**
 * @file   bernoulli_beam_2_example.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Mon Mar 15 2021
 *
 * @brief  Computation of the analytical exemple 1.1 in the TGC vol 6
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "mesh_accessor.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

#define TYPE _bernoulli_beam_2

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
  initialize(argc, argv);
  // Defining the mesh
  Mesh beams(2);

  const auto q = 6000.;
  const auto L = 10.;
  const auto M = -3600.; // Momentum at 3

  auto nb_nodes = 3;
  auto nb_element = nb_nodes - 1;

  MeshAccessor mesh_accessor(beams);
  Array<Real> & nodes = mesh_accessor.getNodes();
  nodes.resize(nb_nodes);

  beams.addConnectivityType(_bernoulli_beam_2);
  Array<UInt> & connectivity = mesh_accessor.getConnectivity(_bernoulli_beam_2);
  connectivity.resize(nb_element);

  nodes.zero();
  nodes(1, 0) = 10;
  nodes(2, 0) = 18;

  for (int i = 0; i < nb_element; ++i) {
    connectivity(i, 0) = i;
    connectivity(i, 1) = i + 1;
  }

  mesh_accessor.makeReady();

  // Defining the materials
  StructuralMechanicsModel model(beams);

  StructuralMaterial mat1;
  mat1.E = 3e10;
  mat1.I = 0.0025;
  mat1.A = 0.01;

  model.addMaterial(mat1);

  StructuralMaterial mat2;
  mat2.E = 3e10;
  mat2.I = 0.00128;
  mat2.A = 0.01;

  model.addMaterial(mat2);

  // Defining the forces
  model.initFull();

  auto & forces = model.getExternalForce();
  auto & displacement = model.getDisplacement();
  auto & boundary = model.getBlockedDOFs();
  const auto & N_M = model.getStress(_bernoulli_beam_2);

  auto & element_material = model.getElementMaterial(_bernoulli_beam_2);

  boundary.set(false);
  forces.zero();
  displacement.zero();

  element_material(1) = 1;

  forces(0, 1) = -q * L / 2.;
  forces(0, 2) = -q * L * L / 12.;
  forces(1, 1) = -q * L / 2.;
  forces(1, 2) = q * L * L / 12.;
  forces(2, 2) = M;
  forces(2, 0) = mat2.E * mat2.A / 18;

  // Defining the boundary conditions
  boundary(0, 0) = true;
  boundary(0, 1) = true;
  boundary(0, 2) = true;
  boundary(1, 1) = true;
  boundary(2, 1) = true;

  model.addDumpFieldVector("displacement");
  model.addDumpField("rotation");
  model.addDumpFieldVector("force");
  model.addDumpField("momentum");

  model.solveStep();

  model.assembleResidual();

  // Post-Processing
  std::cout << " d1 = " << displacement(1, 2) << std::endl;
  std::cout << " d2 = " << displacement(2, 2) << std::endl;
  std::cout << " d3 = " << displacement(1, 0) << std::endl;
  std::cout << " M1 = " << N_M(0, 1) << std::endl;
  std::cout << " M2 = " << N_M(2 * (nb_nodes - 2), 1) << std::endl;

  model.dump();

  finalize();
}
