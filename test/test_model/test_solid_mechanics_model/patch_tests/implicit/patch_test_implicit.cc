/**
 * @file   patch_test_implicit.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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
static types::RMatrix prescribed_strain() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  types::RMatrix strain(spatial_dimension, spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      strain(i, j) = alpha[i][j + 1];
    }
  }
  return strain;
}

template<ElementType type, bool is_plane_strain>
static types::RMatrix prescribed_stress() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  types::RMatrix stress(spatial_dimension, spatial_dimension);

  //plane strain in 2d
  types::RMatrix strain(spatial_dimension, spatial_dimension);
  types::RMatrix pstrain; pstrain = prescribed_strain<type, is_plane_strain>();
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
    stress(0, 0) = E * strain(0, 0);
  } else {
    for (UInt i = 0; i < spatial_dimension; ++i)
      for (UInt j = 0; j < spatial_dimension; ++j) {
	stress(i, j) =  (i == j)*lambda*trace + 2*mu*strain(i, j);
      }
  }

  return stress;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(argc, argv);

  UInt dim = ElementClass<TYPE>::getSpatialDimension();
  const ElementType element_type = TYPE;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;

  std::stringstream filename; filename << TYPE << ".msh";
  mesh_io.read(filename.str(), my_mesh);

  //  MeshUtils::purifyMesh(my_mesh);

  UInt nb_nodes = my_mesh.getNbNodes();

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  if(PLANE_STRAIN)
    my_model.initFull("material_check_stress_plane_strain.dat", _static);
  else
    my_model.initFull("material_check_stress_plane_stress.dat", _static);

  const Vector<Real> & coordinates = my_mesh.getNodes();
  Vector<Real> & displacement = my_model.getDisplacement();
  Vector<bool> & boundary = my_model.getBoundary();

  MeshUtils::buildFacets(my_mesh);
  MeshUtils::buildSurfaceID(my_mesh);

  CSR<UInt> surface_nodes;
  MeshUtils::buildNodesPerSurface(my_mesh, surface_nodes);

  for (UInt s = 0; s < surface_nodes.getNbRows(); ++s) {
    CSR<UInt>::iterator snode = surface_nodes.begin(s);
    for(; snode != surface_nodes.end(s); ++snode) {
      UInt n = *snode;
      std::cout << "Node " << n << std::endl;
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
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  akantu::UInt count = 0;
  my_model.assembleStiffnessMatrix();
  my_model.updateResidual();

  my_model.getStiffnessMatrix().saveMatrix("K.mtx");

  while(!my_model.testConvergenceResidual(2e-4) && (count < 100)) {
    std::cout << "Iter : " << ++count << std::endl;
    my_model.solveStatic();
    my_model.updateResidual();
  }

  if(count > 1) {
    std::cerr << "The code did not converge in 1 step !" << std::endl;
    return EXIT_FAILURE;
  }

  my_model.computeStresses();

  /* ------------------------------------------------------------------------ */
  /* Checks                                                                   */
  /* ------------------------------------------------------------------------ */
  UInt nb_quadrature_points = my_model.getFEM().getNbQuadraturePoints(element_type);

  Vector<Real> & stress_vect = const_cast<Vector<Real> &>(my_model.getMaterial(0).getStress(element_type));
  Vector<Real> & strain_vect = const_cast<Vector<Real> &>(my_model.getMaterial(0).getStrain(element_type));

  Vector<Real>::iterator<types::RMatrix> stress_it = stress_vect.begin(dim, dim);
  Vector<Real>::iterator<types::RMatrix> strain_it = strain_vect.begin(dim, dim);

  types::RMatrix presc_stress; presc_stress = prescribed_stress<TYPE, PLANE_STRAIN>();
  types::RMatrix presc_strain; presc_strain = prescribed_strain<TYPE, PLANE_STRAIN>();

  UInt nb_element = my_mesh.getNbElement(TYPE);

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      types::RMatrix & stress = *stress_it;
      types::RMatrix & strain = *strain_it;

      for (UInt i = 0; i < dim; ++i) {
	for (UInt j = 0; j < dim; ++j) {
	  if(!(std::abs(strain(i, j) - presc_strain(i, j)) < 2e-15)) {
	    std::cerr << "strain[" << i << "," << j << "] = " << strain(i, j) << " but should be = " << presc_strain(i, j) << " (-" << std::abs(strain(i, j) - presc_strain(i, j)) << ") [el : " << el<< " - q : " << q << "]" << std::endl;
	    std::cerr << "computed : " << strain << "reference : " << presc_strain << std::endl;
	    return EXIT_FAILURE;
	  }

	  if(!(std::abs(stress(i, j) - presc_stress(i, j)) < 1e-3)) {
	    std::cerr << "stress[" << i << "," << j << "] = " << stress(i, j) << " but should be = " << presc_stress(i, j) << " (-" << std::abs(stress(i, j) - presc_stress(i, j)) << ") [el : " << el<< " - q : " << q << "]" << std::endl;
	    std::cerr << "computed : " << stress << "reference : " << presc_stress << std::endl;
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

      if(!(std::abs(displacement(n,i) - disp) < 2e-15)) {
	std::cerr << "displacement(" << n << ", " << i <<")=" << displacement(n,i) << " should be equal to " << disp <<  std::endl;
	return EXIT_FAILURE;
      }
    }
  }

  finalize();

  return EXIT_SUCCESS;
}

