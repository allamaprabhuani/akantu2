/**
 * @file   test_hertz_2d.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Apr 30 2019
 * @date last modification: Tue Apr 30 2019
 *
 * @brief  Hertz Test for contact mechanics model class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model.hh"
#include "contact_mechanics_model.hh"
#include "coupler_solid_contact.hh"
#include "non_linear_solver.hh"
#include "dumper_text.hh"
#include "dumper_variable.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  
  UInt max_steps = 5000;
  UInt imposing_steps = 1000;
  UInt damping_interval = 10;
  Real damping_ratio = 0.9;

  Real pressure = 0.01;
  
  std::string mesh_file = "hertz_2d.msh";
  std::string material_file = "material_hertz.dat";

  const UInt spatial_dimension = 2;
  
  initialize(material_file, argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read(mesh_file);

  CouplerSolidContact coupler(mesh);

  auto & solid   = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();

  auto && selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
		     "physical_names",solid);

  solid.initFull(  _analysis_method = _explicit_lumped_mass);
  contact.initFull(_analysis_method = _explicit_dynamic_contact);

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top_body");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "top_body");

  coupler.initFull(_analysis_method = _explicit_dynamic_contact);
  
  Vector<Real> weight = {0, pressure};
  solid.applyBC(BC::Neumann::FromSameDim(weight), "bottom");
  
  Real time_step = solid.getStableTimeStep();
  time_step *= 0.8;

  coupler.setTimeStep(time_step);

  coupler.setBaseName("test-hertz-2d");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("grad_u");
  coupler.addDumpField("stress");

  auto & velocity = solid.getVelocity();

  auto & gaps = contact.getGaps();
  auto & areas = contact.getNodalArea();

  for (UInt s : arange(max_steps)) {

    std::cerr << "Step " << s << std::endl;  
    
    coupler.solveStep();       
        	 
    if (s % damping_interval == 0) {
      for (auto & v : make_view(velocity)) 
	v *= damping_ratio;
    }

    if (s % 100 == 0) {
      coupler.dump();
    }
   
  }

  finalize();
  return EXIT_SUCCESS;

}

