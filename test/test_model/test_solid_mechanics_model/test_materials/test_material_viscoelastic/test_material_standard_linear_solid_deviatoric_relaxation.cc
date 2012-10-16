/**
 * @file   test_material_viscoelastic_relaxation.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Vlad Yastrebov <vladislav.yastrebov@epfl.ch> 
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed May 23 15:54:01 2012
 *
 * @brief  test of the viscoelastic material: relaxation
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
#include <limits>
#include <fstream>
#include <sstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::debug::setDebugLevel(akantu::dblWarning);

  // sim data
  Real T = 10.;
  Real eps = 0.001;

  const UInt dim = 2;
  Real sim_time = 25.;
  Real time_factor = 0.1;

  Real tolerance = 1e-7;

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("test_material_standard_linear_solid_deviatoric_relaxation.msh",mesh);
  const ElementType element_type = _quadrangle_4;
  SolidMechanicsModel model(mesh);

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  model.initFull("material_standard_linear_solid_deviatoric_relaxation.dat", _explicit_dynamic);
  std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();

  std::stringstream filename_sstr;
  filename_sstr << "test_material_standard_linear_solid_deviatoric_relaxation.out";
  std::ofstream output_data;
  output_data.open(filename_sstr.str().c_str());
  output_data << "#[1]-time [2]-sigma_analytic [3+]-sigma_measurements" << std::endl;

  Material & mat = model.getMaterial(0);

  const Vector<Real> & stress = mat.getStress(element_type);

  Real Eta  = mat.getParam<Real>("Eta");
  Real EV   = mat.getParam<Real>("Ev");
  Real Einf = mat.getParam<Real>("Einf");
  Real nu   = mat.getParam<Real>("nu");
  Real Ginf = Einf/(2*(1+nu));
  Real G = EV/(2*(1+nu));
  Real G0 = G + Ginf;
  Real gamma = G/G0;
  Real tau = Eta / EV;
  Real gammainf = Ginf/G0;

  UInt nb_nodes = mesh.getNbNodes();
  const Vector<Real> & coordinate = mesh.getNodes();
  Vector<Real> & displacement = model.getDisplacement();

  /// Setting time step
  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  UInt max_steps = sim_time / time_step;
  UInt out_interval = 1;

  Real time = 0.;

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 0; s <= max_steps; ++s) {

    if(s % 1000 == 0)
      std::cerr << "passing step " << s << "/" << max_steps << std::endl;

    time = s * time_step;
    // impose displacement
    Real epsilon = 0.;
    if (time < T) {
      epsilon = eps * time / T;
    }
    else {
      epsilon = eps;
    }
    for (UInt n=0; n<nb_nodes; ++n) {
      displacement(n,0) = epsilon * coordinate(n,1);
      displacement(n,1) = epsilon * coordinate(n,0);
    }

    // compute stress
    model.updateResidual();

    // print output
    if(s % out_interval == 0) {
      // analytical solution
      Real solution= 0.;
      if (time < T) {
	solution = 2 * G0 * eps / T * (gammainf * time + gamma * tau * (1 - exp(-time/tau)));
      }
      else {
	solution = 2 * G0 * eps * (gammainf + gamma * tau / T * (exp((T-time)/tau) - exp(-time/tau)));
      }
      output_data << s*time_step << " " << solution;

      // data output
      Vector<Real>::const_iterator<types::RMatrix> stress_it = stress.begin(dim, dim);
      Vector<Real>::const_iterator<types::RMatrix> stress_end = stress.end(dim, dim);
      for(;stress_it != stress_end; ++stress_it) {
	output_data << " " << (*stress_it)(0,1) << " " << (*stress_it)(1,0);

	// test error
	Real rel_error_1 = std::abs(((*stress_it)(0,1) - solution) / solution);
	Real rel_error_2 = std::abs(((*stress_it)(1,0) - solution) / solution);
	if (rel_error_1 > tolerance || rel_error_2 > tolerance) {
	  std::cerr << "Relative error: " << rel_error_1 << " " << rel_error_2 << std::endl;
	  return EXIT_FAILURE;
	}
      }
      output_data << std::endl;
    }
  }

  finalize();

  std::cout << "Test successful!" << std::endl;
  return EXIT_SUCCESS;
}
