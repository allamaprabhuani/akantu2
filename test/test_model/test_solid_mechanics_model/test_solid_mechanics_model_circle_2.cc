/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
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

class MyStressFunctor : public SolidMechanicsModel::SurfaceLoadFunctor {
public:
  inline void stress(__attribute__ ((unused)) const types::Vector<Real> & position,
		     types::RMatrix & stress,
		     __attribute__ ((unused)) const types::Vector<Real> & normal,
		     __attribute__ ((unused)) Surface surface_id) {
    stress.eye(1000);
  }
};

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  UInt max_steps = 10000;
  Real epot, ekin;

  Mesh mesh(2);
  MeshIOMSH mesh_io;
  mesh_io.read("circle2.msh", mesh);

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull("material.dat");

  Real time_step = model.getStableTimeStep() / 10.;
  model.setTimeStep(time_step);

  std::cout << "-- Time step : " << time_step << " --" << std::endl;

  model.assembleMassLumped();

  FEM & fem_boundary = model.getFEMBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnControlPoints();

  MyStressFunctor func;
  model.computeForcesFromFunction(func, akantu::_bft_stress);

  model.setBaseName("circle2");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("force"       );
  model.addDumpFieldVector("residual"    );
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.dump();

  for(UInt s = 0; s < max_steps; ++s) {
    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();

    std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;

    if(s % 100 == 0) model.dump();
  }

  return EXIT_SUCCESS;
}



