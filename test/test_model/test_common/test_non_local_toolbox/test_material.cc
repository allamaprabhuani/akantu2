/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 * 
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 * 
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "test_material.hh"

/* -------------------------------------------------------------------------- */
template <Int dim>
TestMaterial<dim>::TestMaterial(SolidMechanicsModel & model, const ID & id)
    : Parent(model, id), grad_u_nl("grad_u non local", *this) {
  this->is_non_local = true;
  this->grad_u_nl.initialize(dim * dim);
}

/* -------------------------------------------------------------------------- */
template <Int dim> void TestMaterial<dim>::registerNonLocalVariables() {
  this->model.getNonLocalManager().registerNonLocalVariable(
      this->gradu.getName(), grad_u_nl.getName(), dim * dim);

  this->model.getNonLocalManager()
      .getNeighborhood(this->getNeighborhoodName())
      .registerNonLocalVariable(grad_u_nl.getName());
}

/* -------------------------------------------------------------------------- */
// Instantiate the material for the 3 dimensions
template class TestMaterial<1>;
template class TestMaterial<2>;
template class TestMaterial<3>;

static bool material_is_allocated_test_material =
    instantiateMaterial<TestMaterial>("test_material");

/* -------------------------------------------------------------------------- */
