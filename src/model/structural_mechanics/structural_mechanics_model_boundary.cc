/**
 * @file   structural_mechanics_model_boundary.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed May 25 15:21:50 2011
 *
 * @brief StructuralMechanicsModel functions to set boundary conditions
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
#include "model.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/**
 * @param myf pointer  to a function that fills a  vector/tensor with respect to
 * passed coordinates
 */
void StructuralMechanicsModel::computeForcesFromFunction(BoundaryFunction myf,
							 BoundaryFunctionType function_type){
  /** function type is
   ** _bft_forces : linear load is given
   ** _bft_stress : stress function is given -> Not already done for this kind of model
   */

  std::stringstream name;
  name << id << ":structuralmechanics:imposed_linear_load";
  Vector<Real> lin_load(0, nb_degree_of_freedom,name.str());
  name.clear();
  name << id << ":structuralmechanics:imposed_stresses";
  Vector<Real> stress_funct(0, nb_degree_of_freedom*nb_degree_of_freedom,name.str());

  UInt offset = 0;
  switch(function_type) {
  case _bft_stress:
    offset = nb_degree_of_freedom * nb_degree_of_freedom; break;
  case _bft_traction:
    offset = nb_degree_of_freedom; break;
  }

  //prepare the loop over element types
  const ElementType type = _bernoulli_beam_2;
  UInt nb_quad           = getFEM().getNbQuadraturePoints(type);
  UInt nb_element        = getFEM().getMesh().getNbElement(type);

  name.clear();
  name << id << ":structuralmechanics:quad_coords";
  Vector<Real> quad_coords(nb_element * nb_quad, spatial_dimension, "quad_coords");


  getFEMClass<MyFEMType>().getShapeFunctions().interpolateOnControlPoints<type>(getFEM().getMesh().getNodes(),
										quad_coords,
										spatial_dimension);
  getFEMClass<MyFEMType>().getShapeFunctions().interpolateOnControlPoints<type>(getFEM().getMesh().getNodes(),
										quad_coords,
										spatial_dimension,
										_not_ghost,
										NULL,
										true,
										0,
										1,
										1);

  Real * imposed_val = NULL;
  switch(function_type) {
  case _bft_stress:
    stress_funct.resize(nb_element*nb_quad);
    imposed_val = stress_funct.values;
    break;
  case _bft_traction:
    lin_load.resize(nb_element*nb_quad);
    imposed_val = lin_load.values;
    break;
  }

  /// sigma/load on each quadrature points
  Real * qcoord = quad_coords.values;
  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quad; ++q) {
      myf(qcoord, imposed_val, NULL, 0);
      imposed_val += offset;
      qcoord += spatial_dimension;
    }
  }

  switch(function_type) {
  case _bft_stress:
    computeForcesByStressTensor(stress_funct,(type)); break;
  case _bft_traction:
    computeForcesByTractionVector(lin_load,(type)); break;
  }
}


/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeForcesByTractionVector(const Vector<Real> & tractions,
							     const ElementType & type){
  AKANTU_DEBUG_IN();
  MyFEMType & fem = getFEMClass<MyFEMType>();
  UInt nb_element = getFEM().getMesh().getNbElement(type);
  UInt nb_nodes_per_element = getFEM().getMesh().getNbNodesPerElement(type);
  UInt nb_quad = getFEM().getNbQuadraturePoints(type);

  // check dimension match
  AKANTU_DEBUG_ASSERT(Mesh::getSpatialDimension(type) == getFEM().getElementDimension(),
		      "element type dimension does not match the dimension of boundaries : " <<
		      getFEM().getElementDimension() << " != " <<
		      Mesh::getSpatialDimension(type));

  // check size of the vector
  AKANTU_DEBUG_ASSERT(tractions.getSize() == nb_quad*nb_element,
		      "the size of the vector should be the total number of quadrature points");

  // check number of components
  AKANTU_DEBUG_ASSERT(tractions.getNbComponent() == nb_degree_of_freedom,
		      "the number of components should be the spatial dimension of the problem");


  Vector<Real> funct(nb_element * nb_quad, nb_degree_of_freedom * nb_nodes_per_element);

  const Vector<Real> & N0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 0);
  const Vector<Real> & M0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 1);
  const Vector<Real> & L0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 2);
  const Vector<Real> & Mp = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 3);
  const Vector<Real> & Lp = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 4);


  //  Vector<Real> n (nb_degree_of_freedom * nb_nodes_per_element, nb_degree_of_freedom);
  funct.clear();

  Real * N0_val = N0.values;
  Real * M0_val = M0.values;
  Real * L0_val = L0.values;
  Real * Mp_val = Mp.values;
  Real * Lp_val = Lp.values;

  types::RMatrix N(nb_degree_of_freedom , nb_degree_of_freedom * nb_nodes_per_element);
  Vector<Real>::iterator< types::RVector> Nt_T = funct.begin(nb_degree_of_freedom * nb_nodes_per_element);
  Vector<Real>::iterator<types::RVector> T = const_cast< Vector<Real> &>(tractions).begin(nb_degree_of_freedom);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quad; ++q) {
      N.clear();
      N(0,0) = N0_val[0];
      N(0,3) = N0_val[1];

      N(1,1) = M0_val[0];
      N(1,2) = L0_val[0];
      N(1,4) = M0_val[1];
      N(1,5) = L0_val[1];

      N(2,1) = Mp_val[0];
      N(2,2) = Lp_val[0];
      N(2,4) = Mp_val[1];
      N(2,5) = Lp_val[1];

      Nt_T->mul<true>(N, *T);

      ++Nt_T;

      N0_val += nb_nodes_per_element;
      M0_val += nb_nodes_per_element;
      L0_val += nb_nodes_per_element;
      Mp_val += nb_nodes_per_element;
      Lp_val += nb_nodes_per_element;
      ++T;
    }

  }

  // allocate the vector that will contain the integrated values
  std::stringstream name;
  name << id << ":solidmechanics:" << type << ":integral_boundary";
  Vector<Real> int_funct(nb_element, nb_degree_of_freedom*nb_nodes_per_element,name.str());
  //do the integration
  getFEM().integrate(funct, int_funct, nb_degree_of_freedom*nb_nodes_per_element, type);
  // assemble the result into force vector
  getFEM().assembleVector(int_funct,*force_momentum,
			  dof_synchronizer->getLocalDOFEquationNumbers(),
			  nb_degree_of_freedom, type);
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeForcesByStressTensor(__attribute__ ((unused)) const Vector<Real> & stresses,
							   __attribute__ ((unused)) const ElementType & type){
  AKANTU_DEBUG_IN();
  /**
  UInt nb_element = getFEMBoundary().getMesh().getNbElement(type);
  UInt nb_quad = getFEMBoundary().getNbQuadraturePoints(type);

  // check dimension match
  AKANTU_DEBUG_ASSERT(Mesh::getSpatialDimension(type) == getFEMBoundary().getElementDimension(),
		      "element type dimension does not match the dimension of boundaries : " <<
		      getFEMBoundary().getElementDimension() << " != " <<
		      Mesh::getSpatialDimension(type));

  // check size of the vector
  AKANTU_DEBUG_ASSERT(stresses.getSize() == nb_quad*nb_element,
		      "the size of the vector should be the total number of quadrature points");

  // check number of components
  AKANTU_DEBUG_ASSERT(stresses.getNbComponent() == spatial_dimension*spatial_dimension,
		      "the number of components should be the dimension of 2-tensors");



  std::stringstream name;
  name << id << ":solidmechanics:" << type << ":traction_boundary";
  Vector<Real> funct(nb_element*nb_quad, spatial_dimension,name.str());
  const Vector<Real> & normals_on_quad = getFEMBoundary().getNormalsOnQuadPoints(type);

  Math::matrix_vector(spatial_dimension,spatial_dimension,stresses,normals_on_quad,funct);
  computeForcesByTractionVector(funct,type);
  */

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
