/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

class SineBoundary : public BC::Dirichlet::DirichletFunctor {
public:
  SineBoundary(Real amp, Real phase, BC::Axis ax = _x)
      : DirichletFunctor(ax), amplitude(amp), phase(phase) {}

public:
  inline void operator()(__attribute__((unused)) UInt node,
                         Vector<bool> & flags, Vector<Real> & primal,
                         const Vector<Real> & coord) const {
    DIRICHLET_SANITY_CHECK;
    flags(axis) = true;
    primal(axis) = -amplitude * std::sin(phase * coord(1));
  }

protected:
  Real amplitude;
  Real phase;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  Int spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("fine_mesh.msh");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  /// boundary conditions
  Vector<Real, 2> traction{.2, .2};
  model.applyBC(SineBoundary(.2, 10., _x), "Fixed_x");
  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "Fixed_y");
  model.applyBC(BC::Neumann::FromTraction(traction), "Traction");

  // output a paraview file with the boundary conditions
  model.setBaseName("plate");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("external_force");
  model.addDumpField("blocked_dofs");
  model.dump();

  return 0;
}
