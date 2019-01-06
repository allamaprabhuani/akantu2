/**
 * @file   tets_phase_field_2d.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Oct 1 2018
 *
 * @brief  test of the class PhaseFieldModel on the 2d square
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "aka_common.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include "phase_field_model.hh"
#include "material.hh"
#include "material_phasefield.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;

/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel &);
void computeStrainOnQuadPoints(SolidMechanicsModel &, PhaseFieldModel &, const GhostType &);
void computeDamageOnQuadPoints(SolidMechanicsModel &, PhaseFieldModel &, const GhostType &);
void gradUToEpsilon(const Matrix<Real> &, Matrix<Real> &);
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {

  std::ofstream os("data.csv");
  os << "#disp stress damage" << std::endl;
  
  initialize("material_coupling.dat", argc, argv);
  
  Mesh mesh(spatial_dimension);
  mesh.read("test_one_element.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

  PhaseFieldModel phase(mesh);
  phase.initFull(_analysis_method = _static);
  
  model.setBaseName("phase_solid");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("damage");
  model.addDumpField("displacement");
  model.dump();
  
  UInt nbSteps = 1000;
  Real increment = 1e-4;

  auto & damage = model.getMaterial(0).getArray<Real>("damage", _quadrangle_4);
  auto & stress = model.getMaterial(0).getArray<Real>("stress", _quadrangle_4);
  
  for (UInt s = 1; s < nbSteps; ++s) {
    applyDisplacement(model);

    model.solveStep();
    computeStrainOnQuadPoints(model, phase, _not_ghost);

    phase.solveStep();
    computeDamageOnQuadPoints(model, phase, _not_ghost);

    os << s*increment << " " << stress(0, 0) << " " << damage(0) << std::endl;
    
    model.dump();
  }

  os.close();
  finalize();

  return EXIT_SUCCESS;

}


/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel & model) {
  auto & displacement = model.getDisplacement();

  auto & positions = model.getMesh().getNodes();
  auto & blocked_dofs = model.getBlockedDOFs();

  
  for (UInt n = 0; n < model.getMesh().getNbNodes(); ++n) {
    if (positions(n, 1) == -0.5) {
      displacement(n, 0) = 0;
      displacement(n, 1) = 0;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n ,1) = true;
    }
    else {
      displacement(n, 0) = 0;
      displacement(n, 1) += 1.e-4;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n ,1) = true;
    }
  }
}

/* -------------------------------------------------------------------------- */
void computeStrainOnQuadPoints(SolidMechanicsModel & solid, PhaseFieldModel & phase, const GhostType & ghost_type) {

  auto & mesh = solid.getMesh();

  auto & strain_on_qpoints = phase.getStrain();
  auto & gradu_on_qpoints  = solid.getMaterial(0).getGradU();
    
  for (auto & type: mesh.elementTypes(spatial_dimension, ghost_type)) {
    auto & strain_on_qpoints_vect = strain_on_qpoints(type, ghost_type);
    auto & gradu_on_qpoints_vect  = gradu_on_qpoints(type, ghost_type);
    for (auto && values:
	   zip(make_view(strain_on_qpoints_vect, spatial_dimension, spatial_dimension),
	       make_view(gradu_on_qpoints_vect,  spatial_dimension, spatial_dimension))) {
      auto & strain = std::get<0>(values);
      auto & grad_u = std::get<1>(values);
      gradUToEpsilon(grad_u, strain);
    }


  }

}


/* -------------------------------------------------------------------------- */
void computeDamageOnQuadPoints(SolidMechanicsModel & solid, PhaseFieldModel & phase, const GhostType & ghost_type) {
  
  auto & fem  = phase.getFEEngine();
  auto & mesh = phase.getMesh();

  switch (spatial_dimension) {
  case 1: {
    auto & mat = static_cast<MaterialPhaseField<1> &>(solid.getMaterial(0));
    auto & damage = mat.getDamage();
    for (auto & type: mesh.elementTypes(spatial_dimension, ghost_type)) {
      auto & damage_on_qpoints_vect = damage(type, ghost_type);
      fem.interpolateOnIntegrationPoints(phase.getDamage(), damage_on_qpoints_vect,
				       1, type, ghost_type); 
    }
    
    break;
  }
  case 2: {
    auto & mat = static_cast<MaterialPhaseField<2> &>(solid.getMaterial(0));
    auto & damage = mat.getDamage();
  
    for (auto & type: mesh.elementTypes(spatial_dimension, ghost_type)) {
      auto & damage_on_qpoints_vect = damage(type, ghost_type);
      fem.interpolateOnIntegrationPoints(phase.getDamage(), damage_on_qpoints_vect,
				       1, type, ghost_type); 
    }
    break;
  }  
  default:
    auto & mat = static_cast<MaterialPhaseField<3> &>(solid.getMaterial(0));
    break;
  }
 
}


/* -------------------------------------------------------------------------- */
void gradUToEpsilon(const Matrix<Real> & grad_u, Matrix<Real> & epsilon) {
  for (UInt i=0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j)
	epsilon(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
  }
}
