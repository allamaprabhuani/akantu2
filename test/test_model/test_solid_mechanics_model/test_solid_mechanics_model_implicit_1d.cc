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
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[])
{
  initialize(argc, argv);
  UInt spatial_dimension = 1;

  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("segment1.msh", mesh);

  SolidMechanicsModel model(mesh);
  model.initFull("material.dat", _static);

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
  UInt count = 0;
  model.updateResidual();
  model.dump();

  while(!model.testConvergenceResidual(1e-1) && (count <= 10)) {
    std::cout << "Iter : " << ++count << std::endl;
    model.assembleStiffnessMatrix();

    model.solveStatic();

    model.getStiffnessMatrix().saveMatrix("Ktmp.mtx");

    model.updateResidual();
    model.computeStresses();
    model.dump();
  }

  finalize();

  return EXIT_SUCCESS;
}
