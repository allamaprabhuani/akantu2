/**
 * @file   test_material_damage.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date   Mon Jan 13 16:05:48 2014
 *
 * @brief  test for material damage class
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
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "element_class.hh"

using namespace akantu;

Real alpha [3][4] = { { 0.01, 0.02, 0.03, 0.04 },
		      { 0.05, 0.06, 0.07, 0.08 },
		      { 0.09, 0.10, 0.11, 0.12 } };

/* -------------------------------------------------------------------------- */
template<ElementType type, bool is_plane_strain>
static Matrix<Real> prescribed_strain() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> strain(spatial_dimension, spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      strain(i, j) = alpha[i][j + 1];
    }
  }
  return strain;
}

template<ElementType type, bool is_plane_strain>
static Matrix<Real> prescribed_stress(Real prescribed_dam) {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> stress(spatial_dimension, spatial_dimension);

  //plane strain in 2d
  Matrix<Real> strain(spatial_dimension, spatial_dimension);
  Matrix<Real> pstrain; pstrain = prescribed_strain<type, is_plane_strain>();
  Real nu = 0.3;
  Real E  = 2.1e11;
  Real trace = 0;

  /// symetric part of the strain tensor
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      strain(i,j) = 0.5 * (pstrain(i, j) + pstrain(j, i));


  for (UInt i = 0; i < spatial_dimension; ++i) trace += strain(i,i);

  Real lambda   = nu * E / ((1 + nu) * (1 - 2*nu));
  Real mu       = E / (2 * (1 + nu));

  if(!is_plane_strain) {
    std::cout << "toto" << std::endl;
    lambda = nu * E / (1 - nu*nu);
  }

  if(spatial_dimension == 1) {
    stress(0, 0) = (1-prescribed_dam) * E * strain(0, 0);
  } else {
    for (UInt i = 0; i < spatial_dimension; ++i)
      for (UInt j = 0; j < spatial_dimension; ++j) {
	stress(i, j) =  (i == j)*(1-prescribed_dam)*lambda*trace + 2*mu*(1-prescribed_dam)*strain(i, j);
      }
  }

  return stress;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(argc, argv);
  
  Real prescribed_dam = 0.1;

  UInt dim = 3;
  const ElementType element_type = _tetrahedron_4;
  const bool plane_strain = true;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;

  std::stringstream filename; filename << "cube_3d_tet_4.msh";
  mesh_io.read(filename.str(), my_mesh);

  //  MeshUtils::purifyMesh(my_mesh);

  UInt nb_nodes = my_mesh.getNbNodes();

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initFull("material_damage.dat", SolidMechanicsModelOptions(_static));

  const Array<Real> & coordinates = my_mesh.getNodes();
  Array<Real> & displacement = my_model.getDisplacement();
  Array<bool> & boundary = my_model.getBoundary();
  MeshUtils::buildFacets(my_mesh);

  my_mesh.createBoundaryGroupFromGeometry();

  // Loop over (Sub)Boundar(ies)
  for(GroupManager::const_element_group_iterator it(my_mesh.element_group_begin());
      it != my_mesh.element_group_end(); ++it) {
    for(ElementGroup::const_node_iterator nodes_it(it->second->node_begin());
        nodes_it!= it->second->node_end(); ++nodes_it) {
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

  /* ------------------------------------------------------------------------ */
  /* Set damage in each gauss point                                           */ 
  /* ------------------------------------------------------------------------ */
  UInt nb_quadrature_points = my_model.getFEM().getNbQuadraturePoints(element_type);

  Array<Real> & dam_vect = const_cast<Array<Real> &>(my_model.getMaterial(0).getInternal("damage")(element_type));
  Array<Real>::iterator<Real> dam_it = dam_vect.begin();
  
  UInt nb_element = my_mesh.getNbElement(element_type);   

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q= 0; q < nb_quadrature_points; ++q) {
      *dam_it = prescribed_dam;

      ++dam_it;
    }
  }
  /* ------------------------------------------------------------------------ */
  /* Static solve                                                             */
  /* ------------------------------------------------------------------------ */
  my_model.solveStep<_scm_newton_raphson_tangent_modified, _scc_residual>(2e-4, 2);

  /* ------------------------------------------------------------------------ */
  /* Checks                                                                   */
  /* ------------------------------------------------------------------------ */
 

  Array<Real> & stress_vect = const_cast<Array<Real> &>(my_model.getMaterial(0).getStress(element_type));

  Array<Real>::iterator< Matrix<Real> > stress_it = stress_vect.begin(dim, dim);

  Matrix<Real> presc_stress; 
  presc_stress = prescribed_stress<element_type, plane_strain>(prescribed_dam);

  Real stress_tolerance = 1e-13;

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & stress = *stress_it;

      Matrix<Real> diff(dim, dim);

      diff  = stress;
      diff -= presc_stress;
      Real stress_error = diff.norm<L_inf>() / stress.norm<L_inf>();

      if(stress_error > stress_tolerance) {
	std::cerr << "stress error: " << stress_error << " > " << stress_tolerance << std::endl;
	std::cerr << "stress: " << stress << std::endl
		  << "prescribed stress: " << presc_stress << std::endl;
	return EXIT_FAILURE;
      } else {
	std::cerr << "stress error: " << stress_error << " < " << stress_tolerance << std::endl;
      }

      ++stress_it;
    }
  }


  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt i = 0; i < dim; ++i) {
      Real disp = alpha[i][0];
      for (UInt j = 0; j < dim; ++j) {
	disp += alpha[i][j + 1] * coordinates(n, j);
      }

      if(!(std::abs(displacement(n,i) - disp) < 2e-15)) {
	std::cerr << "displacement(" << n << ", " << i <<")=" << displacement(n,i) << " should be equal to " << disp <<  std::endl;
	return EXIT_FAILURE;
      }
    }
  }

  finalize();

  return EXIT_SUCCESS;
}


