/**
 * @file   patch_test_explicit.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date   Thu Feb 17 16:05:48 2011
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
#include <iostream>

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "element_class.hh"

using namespace akantu;

Real alpha [3][4] = { { 0.01, 0.02, 0.03, 0.04 },
                      { 0.05, 0.06, 0.07, 0.08 },
                      { 0.09, 0.10, 0.11, 0.12 } };

//Stiffness tensor, rotated by hand
Real C[3][3][3][3] = {{{{95.9445381952, -1.66317905579, 0.831589527639},
                        {-1.66317905579, 65.2849025944, 0.831589527831},
                        {0.831589527639, 0.831589527831, 62.7901340111}},
                       {{-1.66317905579, 21.1509444945, 0.831589527938},
                        {21.1509444945, -1.66317905535, 0.831589527548},
                        {0.831589527938, 0.831589527548, 3.32635811111}},
                       {{0.831589527639, 0.831589527938, 18.6561759111},
                        {0.831589527938, 0.831589527706, 3.32635811116},
                        {18.6561759111, 3.32635811116, -1.66317905552}}},
                      {{{-1.66317905579, 21.1509444945, 0.831589527938},
                        {21.1509444945, -1.66317905535, 0.831589527548},
                        {0.831589527938, 0.831589527548, 3.32635811111}},
                       {{65.2849025944, -1.66317905535, 0.831589527706},
                        {-1.66317905535, 95.9445381937, 0.831589527858},
                        {0.831589527706, 0.831589527858, 62.7901340112}},
                       {{0.831589527831, 0.831589527548, 3.32635811116},
                        {0.831589527548, 0.831589527858, 18.6561759111},
                        {3.32635811116, 18.6561759111, -1.66317905549}}},
                      {{{0.831589527639, 0.831589527938, 18.6561759111},
                        {0.831589527938, 0.831589527706, 3.32635811116},
                        {18.6561759111, 3.32635811116, -1.66317905552}},
                       {{0.831589527831, 0.831589527548, 3.32635811116},
                        {0.831589527548, 0.831589527858, 18.6561759111},
                        {3.32635811116, 18.6561759111, -1.66317905549}},
                       {{62.7901340111, 3.32635811111, -1.66317905552},
                        {3.32635811111, 62.7901340112, -1.66317905549},
                        {-1.66317905552, -1.66317905549, 98.4393067777}}}};


/* -------------------------------------------------------------------------- */
template<ElementType type>
static Matrix<Real> prescribed_grad_u() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> grad_u(spatial_dimension, spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      grad_u(i, j) = alpha[i][j + 1];
    }
  }
  return grad_u;
}

template<ElementType type>
static Matrix<Real> prescribed_stress() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> stress(spatial_dimension, spatial_dimension);

  //plane strain in 2d
  Matrix<Real> strain(spatial_dimension, spatial_dimension);
  Matrix<Real> pstrain; pstrain = prescribed_grad_u<type>();
  Real trace = 0;

  /// symetric part of the strain tensor
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      strain(i,j) = 0.5 * (pstrain(i, j) + pstrain(j, i));


  for (UInt i = 0; i < spatial_dimension; ++i) trace += strain(i,i);

  for (UInt i = 0 ;  i < spatial_dimension ; ++i) {
    for (UInt j = 0 ;  j < spatial_dimension ; ++j) {
      for (UInt k = 0 ;  k < spatial_dimension ; ++k) {
        for (UInt l = 0 ;  l < spatial_dimension ; ++l) {
          stress(i,j) = C[i][j][k][l]*strain(k,l);
        }
      }
    }
  }
  return stress;
}


#define TYPE _tetrahedron_4
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);
  UInt dim = ElementClass<TYPE>::getSpatialDimension();
  const ElementType element_type = TYPE;

  UInt damping_steps = 600000;
  UInt damping_interval = 50;
  Real damping_ratio = 0.99;

  UInt additional_steps = 20000;
  UInt max_steps = damping_steps + additional_steps;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;

  std::stringstream filename; filename << TYPE << ".msh";
  mesh_io.read(filename.str(), my_mesh);

  UInt nb_nodes = my_mesh.getNbNodes();

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initFull("material_anisotropic.dat", _explicit_lumped_mass);


  std::cout << my_model.getMaterial(0) << std::endl;
  Real time_step = my_model.getStableTimeStep()/5.;
  my_model.setTimeStep(time_step);
  my_model.assembleMassLumped();

  std::cout << "The number of time steps is: " << max_steps << " (" << time_step << "s)" << std::endl;

  // boundary conditions
  const Array<Real> & coordinates = my_mesh.getNodes();
  Array<Real> & displacement = my_model.getDisplacement();
  Array<bool> & boundary = my_model.getBoundary();

  MeshUtils::buildFacets(my_mesh);

  my_mesh.getBoundary().createBoundariesFromGeometry();

  // Loop over (Sub)Boundar(ies)
  const Boundary & boundaries = my_mesh.getBoundary();
  for(Boundary::const_iterator it(boundaries.begin()); it != boundaries.end(); ++it) {
    for(akantu::SubBoundary::nodes_const_iterator nodes_it(it->nodes_begin()); nodes_it!= it->nodes_end(); ++nodes_it) {
      UInt n(*nodes_it);
      std::cout << "Node " << *nodes_it << std::endl;
      for (UInt i = 0; i < dim; ++i) {
        displacement(n, i) = alpha[i][0];
        for (UInt j = 0; j < dim; ++j) {
          displacement(n, i) += alpha[i][j + 1] * coordinates(n, j);
        }
        boundary(n, i) = true;
      }
    }
  }

  // Actually, loop over all nodes, since I wanna test a static solution
  for (UInt n = 0 ;  n < nb_nodes ; ++n) {
    for (UInt i = 0 ;  i < dim ; ++i) {
      displacement(n, i) = alpha[i][0];
      for (UInt j = 0 ;  j < dim ; ++j) {
        displacement(n, i) += alpha[i][j+1] * coordinates(n, j);
      }
    }
  }

  Array<Real> & velocity = my_model.getVelocity();

  std::ofstream energy;
  std::stringstream energy_filename; energy_filename << "energy_" << TYPE << ".csv";
  energy.open(energy_filename.str().c_str());
  energy << "id,time,ekin" << std::endl;
  Real ekin_mean = 0.;

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  UInt s;
  for(s = 1; s <= max_steps; ++s) {
    if(s % 10000 == 0) std::cout << "passing step " << s << "/" << max_steps
                                 << " (" << s*time_step << "s)" <<std::endl;

    // damp velocity in order to find equilibrium
    if((s < damping_steps) && (s % damping_interval == 0)) {
      velocity *= damping_ratio;
    }

    if(s % 1000 == 0) {
      ekin_mean = ekin_mean / 1000.;
      std::cout << "Ekin mean = " << ekin_mean << std::endl;
      if (ekin_mean < 1e-10) break;
      ekin_mean = 0.;
    }


    my_model.explicitPred();

    my_model.updateResidual();
    my_model.updateAcceleration();

    my_model.explicitCorr();

    akantu::Real ekin = my_model.getKineticEnergy(); ekin_mean += ekin;

    if(s % 1000 == 0)
      energy << s << "," << s*time_step  << "," << ekin << std::endl;
  }

  energy.close();

  UInt nb_quadrature_points = my_model.getFEM().getNbQuadraturePoints(TYPE);
  Array<Real> & stress_vect = const_cast<Array<Real> &>(my_model.getMaterial(0).getStress(element_type));
  Array<Real> & strain_vect = const_cast<Array<Real> &>(my_model.getMaterial(0).getStrain(element_type));

  Array<Real>::iterator< Matrix<Real> > stress_it = stress_vect.begin(dim, dim);
  Array<Real>::iterator< Matrix<Real> > strain_it = strain_vect.begin(dim, dim);

  Matrix<Real> presc_stress; presc_stress = prescribed_stress<TYPE>();
  Matrix<Real> presc_strain; presc_strain = prescribed_grad_u<TYPE>();

  UInt nb_element = my_mesh.getNbElement(TYPE);

  Real strain_tolerance = 1e-9;
  Real stress_tolerance = 1e2;
  if(s > max_steps) {
    stress_tolerance = 1e4;
    strain_tolerance = 1e-7;
  }


  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & stress = *stress_it;
      Matrix<Real> & strain = *strain_it;

      for (UInt i = 0; i < dim; ++i) {
        for (UInt j = 0; j < dim; ++j) {
          if(!(std::abs(strain(i, j) - presc_strain(i, j)) < strain_tolerance)) {
            std::cerr << "strain[" << i << "," << j << "] = " << strain(i, j) << " but should be = " << presc_strain(i, j) << " (-" << std::abs(strain(i, j) - presc_strain(i, j)) << ") [el : " << el<< " - q : " << q << "]" << std::endl;
            std::cerr << strain << presc_strain << std::endl;
            return EXIT_FAILURE;
          }

          if(!(std::abs(stress(i, j) - presc_stress(i, j)) < stress_tolerance)) {
            std::cerr << "stress[" << i << "," << j << "] = " << stress(i, j) << " but should be = " << presc_stress(i, j) << " (-" << std::abs(stress(i, j) - presc_stress(i, j)) << ") [el : " << el<< " - q : " << q << "]" << std::endl;
            std::cerr << stress << presc_stress << std::endl;
            return EXIT_FAILURE;
          }
        }
      }

      ++stress_it;
      ++strain_it;
    }
  }


  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt i = 0; i < dim; ++i) {
      Real disp = alpha[i][0];
      for (UInt j = 0; j < dim; ++j) {
        disp += alpha[i][j + 1] * coordinates(n, j);
      }

      if(!(std::abs(displacement(n,i) - disp) < 1e-7)) {
        std::cerr << "displacement(" << n << ", " << i <<")=" << displacement(n,i) << " should be equal to " << disp << "(" << displacement(n,i) - disp << ")" <<  std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  finalize();

  return EXIT_SUCCESS;
}
