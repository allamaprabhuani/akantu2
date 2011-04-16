/**
 * @file   test_check_stress.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Wed Feb 16 13:56:42 2011
 *
 * @brief  patch test for elastic material in solid mechanics model
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

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "element_class.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

template<ElementType type>
types::Matrix prescribed_stress();

template<> types::Matrix prescribed_stress<_triangle_3>() {
  UInt spatial_dimension = ElementClass<_triangle_3>::getSpatialDimension();
  types::Matrix stress(spatial_dimension, spatial_dimension);
  stress(0,0) = 0.;
  stress(0,1) = 0.;
  stress(1,0) = 0.;
  stress(1,1) = -2.30769e9;

  return stress;
}

template<ElementType type>
types::Matrix prescribed_strain();

template<> types::Matrix prescribed_strain<_triangle_3>() {
  UInt spatial_dimension = ElementClass<_triangle_3>::getSpatialDimension();
  types::Matrix strain(spatial_dimension, spatial_dimension);
  strain(0,0) = 0.004285714;
  strain(0,1) = 0.;
  strain(1,0) = 0.;
  strain(1,1) = -0.01;

  return strain;
}


int main(int argc, char *argv[])
{
  UInt dim = 2;
  const ElementType element_type = _triangle_3;

  UInt imposing_steps = 1000;
  Real max_displacement = -0.01;

  UInt damping_steps = 400000;
  UInt damping_interval = 50;
  Real damping_ratio = 0.99;

  UInt additional_steps = 20000;
  UInt max_steps = damping_steps + additional_steps;

  std::cout << "The number of time steps is: " << max_steps << std::endl;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("square_check_stress.msh", my_mesh);

  UInt nb_nodes = my_mesh.getNbNodes();

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initVectors();
  // initialize the vectors
  my_model.getForce().clear();
  my_model.getVelocity().clear();
  my_model.getAcceleration().clear();
  my_model.getDisplacement().clear();

  my_model.initModel();
  my_model.readMaterials("material_check_stress.dat");
  my_model.initMaterials();


  Real time_step = my_model.getStableTimeStep()/10.;
  my_model.setTimeStep(time_step);
  my_model.assembleMassLumped();

  // boundary conditions
  Vector<UInt> top_nodes(0, 1);
  Vector<UInt> base_nodes(0, 1);

  const Vector<Real> & coordinates = my_mesh.getNodes();
  Vector<Real> & displacement = my_model.getDisplacement();
  Vector<bool> & boundary = my_model.getBoundary();

  for(UInt n = 0; n < nb_nodes; ++n) {
    Real y_coord = coordinates(n, 1);
    if (y_coord > 0.99) {
      boundary(n, 1) = true;
      top_nodes.push_back(n);
    }
    else if (y_coord < 0.01) {
      boundary(n, 1) = true;
      base_nodes.push_back(n);
    }
  }

  Vector<Real> & velocity = my_model.getVelocity();

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {
    if(s % 10000 == 0) std::cout << "passing step " << s << "/" << max_steps
				 << " (" << s*time_step << "s)" <<std::endl;

    // impose normal displacement
    if(s <= imposing_steps) {
      Real current_displacement = max_displacement / (static_cast<Real>(imposing_steps)) * s;
      for(UInt n = 0; n < top_nodes.getSize(); ++n) {
	UInt node = top_nodes(n);
	displacement(node, 1) = current_displacement;
      }
    }

    // damp velocity in order to find equilibrium
    if((s < damping_steps) && (s % damping_interval == 0)) {
      velocity *= damping_ratio;
    }

    my_model.explicitPred();
    my_model.updateResidual();
    my_model.updateAcceleration();
    my_model.explicitCorr();
  }

  UInt check_element = 0;
  UInt quadrature_point = 0;
  UInt nb_quadrature_points = ElementClass<_triangle_3>::getNbQuadraturePoint();

  Vector<Real> & stress_vect = const_cast<Vector<Real> &>(my_model.getMaterial(0).getStress(element_type));
  Vector<Real> & strain_vect = const_cast<Vector<Real> &>(my_model.getMaterial(0).getStrain(element_type));

  Vector<Real>::iterator<types::Matrix> stress_it = stress_vect.begin(dim, dim);
  Vector<Real>::iterator<types::Matrix> strain_it = strain_vect.begin(dim, dim);
  stress_it += check_element * nb_quadrature_points + quadrature_point;
  strain_it += check_element * nb_quadrature_points + quadrature_point;

  types::Matrix & stress = *stress_it;
  types::Matrix & strain = *strain_it;
  types::Matrix presc_stress; presc_stress = prescribed_stress<_triangle_3>();
  types::Matrix presc_strain; presc_strain = prescribed_strain<_triangle_3>();

  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j) {
      if(!(std::abs(strain(i, j) - presc_strain(i, j)) < 1e-9)) {
	std::cout << "strain[" << i << "," << j << "] = " << strain(i, j) << " but should be = " << presc_strain(i, j) << " " << std::abs(strain(i, j) - presc_strain(i, j)) << std::endl;
	std::cout << strain << presc_strain << std::endl;
	return EXIT_FAILURE;
      }

      if(!(std::abs(stress(i, j) - presc_stress(i, j)) < 1e4)) {
	std::cout << "stress[" << i << "," << j << "] = " << stress(i, j) << " but should be = " << presc_stress(i, j) << " " << std::abs(stress(i, j) - presc_stress(i, j)) << std::endl;
	std::cout << stress << presc_stress << std::endl;
	return EXIT_FAILURE;
      }
    }
  }

  //  finalize();

  return EXIT_SUCCESS;
}
