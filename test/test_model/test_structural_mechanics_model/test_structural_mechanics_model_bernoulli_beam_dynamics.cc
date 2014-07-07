/**
 * @file   test_dynamic.cc
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @date   Wed Apr 23 12:43:39 2014
 *
 * @brief  
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh_struct.hh"
#include "structural_mechanics_model.hh"
#include "material.hh"

#include <iostream>
using namespace akantu;
/* -------------------------------------------------------------------------- */

#define TYPE _bernoulli_beam_2


static Real analytical_solution(Real time, Real L, Real rho, Real E, Real A, Real I, Real F) {
  Real omega = M_PI*M_PI/L/L*sqrt(E*I/rho);
  Real sum = 0.;
  UInt  i=5;
  for(UInt n=1; n <= i; n +=2) {
    sum += (1. - cos(n*n*omega*time))/pow(n, 4);
  }

  return 2.*F*pow(L, 3)/pow(M_PI, 4)/E/I*sum;
}


//load 
const Real F = 0.5e4;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){
  initialize(argc, argv);
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);
  const ElementType type = _bernoulli_beam_2;

  /* -------------------------------------------------------------------------- */
  // Mesh
  UInt nb_element = 8;
  UInt nb_nodes = nb_element + 1;
  Real total_length = 10.;
  Real length = total_length/nb_element;
  Real heigth = 1.;
  
  Array<Real> & nodes = const_cast<Array<Real> &>(beams.getNodes());
  nodes.resize(nb_nodes);
  
  beams.addConnectivityType(type);
  Array<UInt> &connectivity = const_cast<Array<UInt> &>(beams.getConnectivity(type));
  connectivity.resize(nb_element);

  beams.initNormals();
  Array<Real> & normals = const_cast<Array<Real> &>(beams.getNormals(type));
  normals.resize(nb_element);
  
  for (UInt i = 0; i < nb_nodes; ++i){
    nodes(i,0) = i * length;
    nodes(i,1) = 0;
  }
  
  for (UInt i = 0; i < nb_element; ++i){
    connectivity(i,0) = i;
    connectivity(i,1) = i+1;
  }
  
  /* -------------------------------------------------------------------------- */
  // Materials
  
  StructuralMechanicsModel model(beams);
  
  StructuralMaterial mat1;
  mat1.E = 120e6;
  mat1.rho = 1000;
  mat1.A = heigth;
  mat1.I = heigth*heigth*heigth/12;
  model.addMaterial(mat1);

  /* -------------------------------------------------------------------------- */
  // Forces
  // model.initFull();
  model.initFull(StructuralMechanicsModelOptions(_implicit_dynamic));

  const Array<Real> &position = beams.getNodes();
  Array<Real> &forces = model.getForce();
  Array<Real> &displacement = model.getDisplacement();
  Array<bool> &boundary = model.getBlockedDOFs();

  UInt node_to_print = -1;

  forces.clear();
  displacement.clear();
  //  boundary.clear();
  //model.getElementMaterial(type)(i,0) = 0;
  //model.getElementMaterial(type)(i,0) = 1;
  for (UInt i = 0; i < nb_element; ++i) {
    model.getElementMaterial(type)(i,0) = 0;
  }

  for (UInt n =0; n<nb_nodes; ++n){
    Real x = position(n, 0);
    //    Real y = position(n, 1);

    if (Math::are_float_equal(x, total_length/2.)){
      forces(n,1) = F;
      node_to_print = n;
    }
  }

  std::ofstream pos;
  pos.open("position.csv");
  if(!pos.good()) {
    std::cerr << "Cannot open file" << std::endl;
    exit(EXIT_FAILURE);
  }
  pos << "id,time,position,solution" << std::endl;

  // model.computeForcesFromFunction<type>(load, _bft_traction)
  /* -------------------------------------------------------------------------- */
  // Boundary conditions
  // true ~ displacement is blocked
  boundary(0,0) = true;
  boundary(0,1) = true;
  boundary(nb_nodes-1,1) = true;
  /* -------------------------------------------------------------------------- */
  // "Solve"
  
  Real time = 0;
  model.assembleStiffnessMatrix();
  model.assembleMass();
  model.dump();
  model.getStiffnessMatrix().saveMatrix("K.mtx");
  model.getMassMatrix().saveMatrix("M.mt");
  Real time_step = 1e-4;
  model.setTimeStep(time_step);

  std::cout << "Time" << "  |   " << "Mid-Span Displacement" << std::endl;

  /// time loop
  for (UInt s = 1; time < 0.64; ++s) {

    model.solveStep<_scm_newton_raphson_tangent, _scc_increment>(1e-12, 1000);
  
    pos << s << "," << time << "," << displacement(node_to_print, 1) << "," << analytical_solution(s*time_step, total_length, mat1.rho, mat1.E, mat1.A, mat1.I, F) << std::endl;
    //    pos << s << "," << time << "," << displacement(node_to_print, 1) << "," << analytical_solution(s*time_step) << std::endl;

    time += time_step;
    if(s % 100 == 0) 
      std::cout << time << "  |   " << displacement(node_to_print, 1) << std::endl;  
    model.dump();
  }

  pos.close();

  finalize();

  return EXIT_SUCCESS;
}

