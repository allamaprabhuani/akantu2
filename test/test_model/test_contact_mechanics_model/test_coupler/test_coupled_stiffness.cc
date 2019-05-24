/**
 * @file   test_contact_mechanics_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Apr 30 2019
 * @date last modification: Tue Apr 30 2019
 *
 * @brief  Test for contact mechanics model class
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
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  
  const UInt spatial_dimension = 2;
  initialize("material_stiffness.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("flat_on_flat.msh");
  
  CouplerSolidContact coupler(mesh);

  auto & solid   = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();

  solid.initFull(  _analysis_method = _static);
  contact.initFull(_analysis_method = _implicit_contact);

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top_body");
  solid.applyBC(BC::Dirichlet::IncrementValue(-0.001, _y), "top_body");

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");
  
  coupler.initFull(_analysis_method = _implicit_contact);

  auto & solver = coupler.getNonLinearSolver();
  solver.set("max_iterations", 1);
  solver.set("threshold", 1e-1);
  solver.set("convergence_type", _scc_solution);
  
  coupler.setBaseName("test-coupled-stiffness");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("grad_u");
  coupler.addDumpField("stress");

  
  auto & before_assembly =
    const_cast<SparseMatrix &>(coupler.getDOFManager().getNewMatrix("K", _symmetric));
  
  solid.assembleStiffnessMatrix();

  auto & solid_assembly =
    const_cast<SparseMatrix &>(coupler.getDOFManager().getMatrix("K"));
  solid_assembly.saveMatrix("solid_assembly.mtx");

  auto & displacement = solid.getDisplacement();
   
  contact.search(displacement);
  contact.assembleStiffnessMatrix();

  auto contact_map = contact.getContactMap();
  auto nb_contacts = contact_map.size();
				     
  auto & contact_assembly =
    const_cast<SparseMatrix &>(coupler.getDOFManager().getMatrix("K"));
  contact_assembly.saveMatrix("contact_assembly.mtx");

  solid.assembleInternalForces();
  contact.assembleInternalForces();

  coupler.dump();
  
  Array<Real> & contact_force = contact.getInternalForce();

  for (UInt n : arange(contact_force.size())) {
    std::cerr << contact_force(n, 1) << std::endl;
  }
  
  if (solid_assembly.size() == contact_assembly.size() and nb_contacts > 0) {
    std::cerr << "size of stiffness matrix of solid   = " << solid_assembly.size() << std::endl;
    std::cerr << "size of stiffness matrix of coupled = " << contact_assembly.size() << std::endl;
    std::cerr << "number of contacts = " << nb_contacts << std::endl;

    for (auto & pair : contact_map) {
      std::cerr << "Node " << pair.first << " in contact with " << pair.second.master <<
	" of gap " <<  pair.second.gap  << std::endl;
    }
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}

