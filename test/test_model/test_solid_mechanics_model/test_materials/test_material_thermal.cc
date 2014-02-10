/**
 * @file   test_solid_mechanics_model_square.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Sep 22 13:39:02 2010
 *
 * @brief  test of the class SolidMechanicsModel
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "fem.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[])
{
  debug::setDebugLevel(dblWarning);
  initialize(argc, argv);

  Math::setTolerance(1.e-13);

  Mesh mesh(2);
  MeshIOMSH mesh_io;
  mesh_io.read("square.msh", mesh);

  SolidMechanicsModel model(mesh);
  model.initFull("material_thermal.dat", SolidMechanicsModelOptions(_static));

  mesh.computeBoundingBox();
  Real xmin = mesh.getXMin();
  Real xmax = mesh.getXMax();
  Real ymin = mesh.getYMin();
  Real ymax = mesh.getYMax();

  Array<Real> & pos = mesh.getNodes();
  Array<bool> & boundary = model.getBoundary();
  Array<Real> & disp = model.getDisplacement();

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (Math::are_float_equal(pos(i, 0), xmin)) {
      boundary(i, 0) = true;
    }

    if (Math::are_float_equal(pos(i, 1), ymin)) {
      boundary(i, 1) = true;
    }
  }

  model.setBaseName("test_material_thermal");
  model.addDumpField("displacement");
  model.addDumpField("strain");
  model.addDumpField("stress");
  model.addDumpField("delta_T");

  model.assembleStiffnessMatrix();
  model.updateResidual();

  model.solveStatic();
  model.updateResidual();

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (Math::are_float_equal(pos(i, 0), xmax) && Math::are_float_equal(pos(i, 1), ymax)) {
      if (!Math::are_float_equal(disp(i, 0), 1.0) || !Math::are_float_equal(disp(i, 1), 1.0)) {
	AKANTU_DEBUG_ERROR("Test not passed");
        return EXIT_FAILURE;
      }
    }
  }

  model.dump();

  finalize();

  std::cout << "Test passed" << std::endl;
  return EXIT_SUCCESS;
}



