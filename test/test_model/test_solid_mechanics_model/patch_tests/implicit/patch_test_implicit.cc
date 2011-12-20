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
#  include "io_helper.hh"
using namespace iohelper;
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

Real alpha [3][4] = { { 0.01, 0.02, 0.03, 0.04 },
		      { 0.05, 0.06, 0.07, 0.08 },
		      { 0.09, 0.10, 0.11, 0.12 } };


/* -------------------------------------------------------------------------- */
template<ElementType type>
types::Matrix prescribed_strain() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  types::Matrix strain(spatial_dimension, spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      strain(i, j) = alpha[i][j + 1];
    }
  }
  return strain;
}

template<ElementType type>
types::Matrix prescribed_stress() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  types::Matrix stress(spatial_dimension, spatial_dimension);

  //plane strain in 2d
  types::Matrix strain(spatial_dimension, spatial_dimension);
  types::Matrix pstrain; pstrain = prescribed_strain<type>();
  Real nu = 0.3;
  Real E  = 2.1e11;
  Real trace = 0;

  /// symetric part of the strain tensor
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      strain(i,j) = 0.5 * (pstrain(i, j) + pstrain(j, i));


  for (UInt i = 0; i < spatial_dimension; ++i) trace += strain(i,i);

  Real Ep = E / (1 + nu);
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      stress(i, j) = Ep * strain(i,j);
      if(i == j) stress(i, j) += Ep * (nu / (1 - 2*nu)) * trace;
    }

  return stress;
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
void paraviewInit(Dumper & dumper, const SolidMechanicsModel & model);
void paraviewDump(Dumper & dumper);
#endif

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

  my_model.initImplicit();

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
  my_model.updateResidual();

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  paraviewInit(dumper, my_model);
#endif

  while(!my_model.testConvergenceResidual(2e-4) && (count < 100)) {
    std::cout << "Iter : " << ++count << std::endl;

    my_model.assembleStiffnessMatrix();
    my_model.solveStatic();
    my_model.updateResidual();
  }

#ifdef AKANTU_USE_IOHELPER
  paraviewDump(dumper);
#endif

  if(count > 1) {
    std::cerr << "The code did not converge in 1 step !" << std::endl;
    return EXIT_FAILURE;
  }


  /* ------------------------------------------------------------------------ */
  /* Checks                                                                   */
  /* ------------------------------------------------------------------------ */
  UInt nb_quadrature_points = my_model.getFEM().getNbQuadraturePoints(element_type);

  Vector<Real> & stress_vect = const_cast<Vector<Real> &>(my_model.getMaterial(0).getStress(element_type));
  Vector<Real> & strain_vect = const_cast<Vector<Real> &>(my_model.getMaterial(0).getStrain(element_type));

  Vector<Real>::iterator<types::Matrix> stress_it = stress_vect.begin(dim, dim);
  Vector<Real>::iterator<types::Matrix> strain_it = strain_vect.begin(dim, dim);

  types::Matrix presc_stress; presc_stress = prescribed_stress<TYPE>();
  types::Matrix presc_strain; presc_strain = prescribed_strain<TYPE>();

  UInt nb_element = my_mesh.getNbElement(TYPE);

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      types::Matrix & stress = *stress_it;
      types::Matrix & strain = *strain_it;

      for (UInt i = 0; i < dim; ++i) {
	for (UInt j = 0; j < dim; ++j) {
	  if(!(std::abs(strain(i, j) - presc_strain(i, j)) < 1e-15)) {
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

      if(!(std::abs(displacement(n,i) - disp) < 1e-15)) {
	std::cerr << "displacement(" << n << ", " << i <<")=" << displacement(n,i) << " should be equal to " << disp <<  std::endl;
	return EXIT_FAILURE;
      }
    }
  }

  // std::cout << "Strain : " << strain;
  // std::cout << "Stress : " << stress;

  //  finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/* Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
template <ElementType type> ElemType paraviewType();

template <> ElemType paraviewType<_segment_2>()      { return LINE1; };
template <> ElemType paraviewType<_segment_3>()      { return LINE2; };
template <> ElemType paraviewType<_triangle_3>()     { return TRIANGLE1; };
template <> ElemType paraviewType<_triangle_6>()     { return TRIANGLE2; };
template <> ElemType paraviewType<_quadrangle_4>()   { return QUAD1; };
template <> ElemType paraviewType<_quadrangle_8>()   { return QUAD2; };
template <> ElemType paraviewType<_tetrahedron_4>()  { return TETRA1; };
template <> ElemType paraviewType<_tetrahedron_10>() { return TETRA2; };
template <> ElemType paraviewType<_hexahedron_8>()   { return HEX1; };

void paraviewInit(Dumper & dumper, const SolidMechanicsModel & model) {
  UInt spatial_dimension = ElementClass<TYPE>::getSpatialDimension();
  UInt nb_nodes   = model.getFEM().getMesh().getNbNodes();
  UInt nb_element = model.getFEM().getMesh().getNbElement(TYPE);

  std::stringstream filename; filename << "out_" << TYPE;

  dumper.SetMode(TEXT);
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, filename.str().c_str());
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(TYPE).values,
			 paraviewType<TYPE>(), nb_element, C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model.getMass().values,
			  1, "mass");
  dumper.AddNodeDataField(model.getForce().values,
			  spatial_dimension, "applied_force");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(TYPE).values,
   			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(TYPE).values,
   			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(Dumper & dumper) {
  dumper.Dump();
}

#endif
