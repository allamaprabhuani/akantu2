/**
 * @file   patch_test_linear_elastic_implicit.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jan 30 2018
 *
 * @brief  Patch test for SolidMechanics implicit
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "patch_test_linear_solid_mechanics_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

TYPED_TEST(TestPatchTestSMMLinear, Implicit) {
  std::string filename = "material_check_stress_plane_stress.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain.dat";

  this->initModel(_implicit_dynamic, filename);

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();
  // set the position of all nodes to the static solution
  for (auto && tuple : zip(make_view(coordinates, this->dim),
                           make_view(displacement, this->dim))) {
    this->setLinearDOF(std::get<1>(tuple), std::get<0>(tuple));
  }
  for (UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, 1e-16);

  this->checkAll();
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestSMMLinear, Static) {
  std::string filename = "material_check_stress_plane_stress.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain.dat";

  this->initModel(_static, filename);

  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  this->model->solveStep();

  this->checkAll();
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestSMMLinear, ImplicitFiniteDeformation) {
  std::string filename = "material_check_stress_plane_stress_finite_deformation.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain_finite_deformation.dat";

  this->initModel(_implicit_dynamic, filename);

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();
  // set the position of all nodes to the static solution
  for (auto && tuple : zip(make_view(coordinates, this->dim),
                           make_view(displacement, this->dim))) {
    this->setLinearDOF(std::get<1>(tuple), std::get<0>(tuple));
  }
  for (UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, 1e-16);

  this->checkAll();
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestSMMLinear, StaticFiniteDeformation) {
  std::string filename = "material_check_stress_plane_stress_finite_deformation.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain_finite_deformation.dat";

  this->initModel(_static, filename);

  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  this->model->solveStep();

  this->checkAll();
}
