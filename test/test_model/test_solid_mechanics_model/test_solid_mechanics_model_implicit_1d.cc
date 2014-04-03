/**
 * @file   test_solid_mechanics_model_implicit_1d.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Mar 03 16:09:49 2011
 *
 * @brief  test of traction in implicit
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
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[])
{
  initialize("material.dat", argc, argv);
  UInt spatial_dimension = 1;

  Mesh mesh(spatial_dimension);
  mesh.read("segment1.msh");

  SolidMechanicsModel model(mesh);
  model.initFull( SolidMechanicsModelOptions(_static));

  std::cout << model.getMaterial(0) << std::endl;

  /// boundary conditions
  model.getBoundary()(0,0) = true;
  model.getForce()(1,0) = 1000;

  model.setBaseName("implicit_1d");
  model.addDumpField("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );

  debug::setDebugLevel(dblInfo);

  model.dump();

  model.solveStep<_scm_newton_raphson_tangent_modified, _scc_residual>(1e-1, 10);

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
