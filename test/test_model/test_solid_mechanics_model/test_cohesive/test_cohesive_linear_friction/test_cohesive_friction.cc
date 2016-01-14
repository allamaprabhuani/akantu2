/**
 * @file   test_cohesive_friction.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date   Wed Jan 13 11:56:18 2016
 *
 * @brief testing the correct behavior of the friction law included in
 * the cohesive linear law, in implicit
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or  modify it under the
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
#include "mesh_io.hh"
#include "mesh_utils.hh"
#include "static_communicator.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material.hh"
#include "material_cohesive.hh"
#include "aka_grid_dynamic.hh"
#include "dof_synchronizer.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>

#include <iomanip>

/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {

  initialize("material.dat", argc, argv);

  Math::setTolerance(1.e-15);
  UInt spatial_dimension = 2;
  const ElementType type = _cohesive_2d_4;

  Mesh mesh(spatial_dimension);

  mesh.read("mesh_cohesive_friction.msh");
  
  // Create the model
  SolidMechanicsModelCohesive model(mesh);

  // Model initialization
  model.initFull(SolidMechanicsModelCohesiveOptions(_static, true));

  //  CohesiveElementInserter inserter(mesh);
  model.limitInsertion(_y, -0.001, 0.001);
  model.updateAutomaticInsertion();

  Real eps = 1e-10;
  Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();
  const Array<Real> & residual = model.getResidual();
  Array<Real> & cohe_opening = const_cast<Array<Real> &>(model.getMaterial("interface").getInternal<Real>("opening")(type));
  Array<Real> & friction_force = const_cast<Array<Real> &>(model.getMaterial("interface").getInternal<Real>("friction_force")(type));

  // Boundary conditions
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if(pos(i,1) < -0.49 || pos(i,1) > 0.49){
      boun(i,0) = true;
      boun(i,1) = true;
    }
  }

  bool passed = true;
  Real tolerance = 1e-13;
  Real error;
  bool load_reduction = false;
  Real tol_increase_factor = 1e5;
  Real increment = 1.0e-4;

  model.synchronizeBoundaries();
  model.updateResidual();

  /* -------------------------------------------- */
  /* LOADING PHASE to introduce cohesive elements */
  /* -------------------------------------------- */

  for (UInt nstep = 0; nstep < 100; ++nstep){

    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (pos(n,1) > 0.49)
	disp(n,1) += increment;
    }

    model.solveStepCohesive<_scm_newton_raphson_tangent, _scc_increment>(tolerance, error, 25, load_reduction, tol_increase_factor);

    if (error > tolerance){
      AKANTU_DEBUG_ERROR("Convergence not reached in the mode I loading phase");
      passed = false;
    }

  }

  /* --------------------------------------------------------- */
  /* UNLOADING PHASE to bring cohesive elements in compression */
  /* --------------------------------------------------------- */

  for (UInt nstep = 0; nstep < 110; ++nstep){

    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (pos(n,1) > 0.49)
	disp(n,1) -= increment;
    }

    model.solveStepCohesive<_scm_newton_raphson_tangent, _scc_increment>(tolerance, error, 25, load_reduction, tol_increase_factor);

    if (error > tolerance){
      AKANTU_DEBUG_ERROR("Convergence not reached in the mode I unloading phase");
      passed = false;
    }

  }


  /* -------------------------------------------------- */
  /* SHEAR PHASE - displacement towards right           */
  /* -------------------------------------------------- */

  increment *= 2;

  for (UInt nstep = 0; nstep < 30; ++nstep){

    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (pos(n,1) > 0.49)
	disp(n,0) += increment;
    }

    model.solveStepCohesive<_scm_newton_raphson_tangent, _scc_increment>(tolerance, error, 25, load_reduction, tol_increase_factor);

    if (error > tolerance){
      AKANTU_DEBUG_ERROR("Convergence not reached in the shear loading phase");
      passed = false;
    }

  }

  /* ---------------------------------------------------*/
  /* Check the horizontal component of the reaction     */
  /* ---------------------------------------------------*/

  // Friction + mode II cohesive behavior
  Real reac_X = 0.;

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (pos(i,1) > 0.49)
      reac_X += residual(i,0);
  }
  if (std::abs(reac_X - (-23.093432803069547)) > eps)
    passed = false;

  // Only friction
  Real friction = friction_force(0,0) + friction_force(1,0);

  if (std::abs(friction - (-21.714260649226127)) > eps)
    passed = false;


  /* -------------------------------------------------- */
  /* SHEAR PHASE - displacement back to zero            */
  /* -------------------------------------------------- */

  for (UInt nstep = 0; nstep < 30; ++nstep){

    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (pos(n,1) > 0.49)
	disp(n,0) -= increment;
    }

    model.solveStepCohesive<_scm_newton_raphson_tangent, _scc_increment>(tolerance, error, 25, load_reduction, tol_increase_factor);

    if (error > tolerance){
      AKANTU_DEBUG_ERROR("Convergence not reached in the shear unloading phase");
      passed = false;
    }

  }


  /* ------------------------------------------------------- */
  /* Check the horizontal component of the reaction and      */
  /* the residual relative sliding in the cohesive elements  */
  /* ------------------------------------------------------- */

  // Friction + mode II cohesive behavior
  reac_X = 0.;

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (pos(i,1) > 0.49)
      reac_X += residual(i,0);
  }
  if (std::abs(reac_X - 21.517959089220852) > eps)
    passed = false;

  // Only friction
  friction = 0.;
  friction = friction_force(0,0) + friction_force(1,0);

  if (std::abs(friction - 21.731377046498903) > eps)
    passed = false;

  // Residual sliding
  Real sliding[2];
  sliding[0] = cohe_opening(0,0);
  sliding[1] = cohe_opening(1,0);

  if (std::abs(sliding[0] - (-0.0007678018950562827)) > eps)
    passed = false;


  if (passed)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;

  finalize();

  return EXIT_SUCCESS;
}
