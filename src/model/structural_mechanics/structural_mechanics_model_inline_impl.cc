/**
 * @file   structural_mechanics_model_inline_impl.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu May  5 19:48:07 2011
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
template<ElementType type>
inline UInt StructuralMechanicsModel::getTangentStiffnessVoigtSize() {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return 0;
}

template<>
inline UInt StructuralMechanicsModel::getTangentStiffnessVoigtSize<_bernoulli_beam_2>() {
  return 2;
}

template<>
inline UInt StructuralMechanicsModel::getTangentStiffnessVoigtSize<_bernoulli_beam_3>() {
  return 4;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  SparseMatrix & K = *stiffness_matrix;

  UInt nb_element                 = getFEM().getMesh().getNbElement(type);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(type);

  UInt tangent_size = getTangentStiffnessVoigtSize<type>();

  Vector<Real> * tangent_moduli =
    new Vector<Real>(nb_element * nb_quadrature_points, tangent_size * tangent_size,
		     "tangent_stiffness_matrix");

  tangent_moduli->clear();

  computeTangentModuli<type>(*tangent_moduli);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = nb_degree_of_freedom * nb_nodes_per_element;

  Vector<Real> * bt_d_b = new Vector<Real>(nb_element*nb_quadrature_points,
					   bt_d_b_size * bt_d_b_size,
					   "B^t*D*B");

  Vector<Real> * b = new Vector<Real>(nb_element*nb_quadrature_points,
				      tangent_size*bt_d_b_size,
				      "B");

  transferBMatrixToSymVoigtBMatrix<type>(*b);

  types::RMatrix Bt_D(bt_d_b_size, tangent_size);
  types::RMatrix BT(tangent_size, bt_d_b_size);

  Vector<Real>::iterator<types::RMatrix> B = b->begin(tangent_size, bt_d_b_size);
  Vector<Real>::iterator<types::RMatrix> D = tangent_moduli->begin(tangent_size, tangent_size);
  Vector<Real>::iterator<types::RMatrix> Bt_D_B = bt_d_b->begin(bt_d_b_size, bt_d_b_size);
  Vector<Real>::iterator<types::RMatrix> T = rotation_matrix(type).begin(bt_d_b_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e, ++T) {
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++B, ++D, ++Bt_D_B) {
      BT.mul<false, false>(*B, *T);
      Bt_D.mul<true, false>(BT, *D);
      Bt_D_B->mul<false, false>(Bt_D, BT);
    }
  }

  delete b;
  delete tangent_moduli;

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Vector<Real> * int_bt_d_b = new Vector<Real>(nb_element,
					   bt_d_b_size * bt_d_b_size,
					   "int_B^t*D*B");

  getFEM().integrate(*bt_d_b, *int_bt_d_b,
		     bt_d_b_size * bt_d_b_size,
		     type);

  delete bt_d_b;

  getFEM().assembleMatrix(*int_bt_d_b, K, nb_degree_of_freedom, type);

  delete int_bt_d_b;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
void StructuralMechanicsModel::computeTangentModuli(Vector<Real> & tangent_moduli) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


/* -------------------------------------------------------------------------- */
template<ElementType type>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix(Vector<Real> & b, bool local) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
void StructuralMechanicsModel::computeStressOnQuad() {
  AKANTU_DEBUG_IN();

  Vector<Real> & sigma  = stress(type, _not_ghost);

  sigma.clear();
  const Mesh & mesh = getFEM().getMesh();

  UInt nb_element                 = mesh.getNbElement(type);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(type);

  UInt tangent_size = getTangentStiffnessVoigtSize<type>();

  Vector<Real> * tangent_moduli =
    new Vector<Real>(nb_element*nb_quadrature_points, tangent_size * tangent_size,
		     "tangent_stiffness_matrix");

  computeTangentModuli<type>(*tangent_moduli);

  /// compute DB
  UInt d_b_size = nb_degree_of_freedom * nb_nodes_per_element;

  Vector<Real> * d_b = new Vector<Real>(nb_element*nb_quadrature_points,
					d_b_size * tangent_size,
					"D*B");

  Vector<Real> * b = new Vector<Real>(nb_element*nb_quadrature_points,
				      tangent_size*d_b_size,
				      "B");

  transferBMatrixToSymVoigtBMatrix<type>(*b);

  Vector<Real>::iterator<types::RMatrix> B = b->begin(tangent_size, d_b_size);
  Vector<Real>::iterator<types::RMatrix> D = tangent_moduli->begin(tangent_size, tangent_size);
  Vector<Real>::iterator<types::RMatrix> D_B = d_b->begin(tangent_size, d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++B, ++D, ++D_B) {
      D_B->mul<false, false>(*D, *B);
    }
  }

  delete b;
  delete tangent_moduli;

  /// compute DBu
  D_B = d_b->begin(tangent_size, d_b_size);
  Vector<Real>::iterator< types::RVector> DBu = sigma.begin(tangent_size);
  types::RVector ul (d_b_size);

  Vector<Real> u_el(0, d_b_size);
  FEM::extractNodalToElementField(mesh, *displacement_rotation, u_el, type);

  Vector<Real>::iterator<types::RVector> ug = u_el.begin(d_b_size);
  Vector<Real>::iterator<types::RMatrix> T = rotation_matrix(type).begin(d_b_size, d_b_size);

  for (UInt e = 0; e < nb_element; ++e, ++T, ++ug) {
    ul.mul<false>(*T, *ug);
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_B, ++DBu) {
      DBu->mul<false>(*D_B, ul);
    }
  }

  delete d_b;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
void StructuralMechanicsModel::computeForcesByLocalTractionVector(const Vector<Real> & tractions) {
  AKANTU_DEBUG_IN();

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


  Vector<Real> Nvoigt(nb_element * nb_quad, nb_degree_of_freedom * nb_degree_of_freedom * nb_nodes_per_element);
  transferNMatrixToSymVoigtNMatrix<type>(Nvoigt);

  Vector<Real>::const_iterator<types::RMatrix> N_it = Nvoigt.begin(nb_degree_of_freedom,
								   nb_degree_of_freedom * nb_nodes_per_element);
  Vector<Real>::const_iterator<types::RMatrix> T_it = rotation_matrix(type).begin(nb_degree_of_freedom * nb_nodes_per_element,
										  nb_degree_of_freedom * nb_nodes_per_element);
  Vector<Real>::const_iterator<types::RVector> te_it = tractions.begin(nb_degree_of_freedom);

  Vector<Real> funct(nb_element * nb_quad, nb_degree_of_freedom * nb_nodes_per_element, 0.);
  Vector<Real>::iterator< types::RVector> Fe_it = funct.begin(nb_degree_of_freedom * nb_nodes_per_element);

  types::Vector<Real> fe(nb_degree_of_freedom * nb_nodes_per_element);
  for (UInt e = 0; e < nb_element; ++e, ++T_it) {
    const types::RMatrix & T = *T_it;
    for (UInt q = 0; q < nb_quad; ++q, ++N_it, ++te_it, ++Fe_it) {
      const types::RMatrix & N = *N_it;
      const types::RVector & te = *te_it;
      types::RVector & Fe = *Fe_it;

      // compute N^t tl
      fe.mul<true>(N, te);
      // turn N^t tl back in the global referential
      Fe.mul<true>(T, fe);
    }
  }

  // allocate the vector that will contain the integrated values
  std::stringstream name;
  name << id << type << ":integral_boundary";
  Vector<Real> int_funct(nb_element, nb_degree_of_freedom * nb_nodes_per_element, name.str());

  //do the integration
  getFEM().integrate(funct, int_funct, nb_degree_of_freedom*nb_nodes_per_element, type);

  // assemble the result into force vector
  getFEM().assembleVector(int_funct,*force_momentum,
			  dof_synchronizer->getLocalDOFEquationNumbers(),
			  nb_degree_of_freedom, type);
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<ElementType type>
void StructuralMechanicsModel::computeForcesByGlobalTractionVector(const Vector<Real> & traction_global){
 AKANTU_DEBUG_IN();
  UInt nb_element = getFEM().getMesh().getNbElement(type);
  UInt nb_quad = getFEM().getNbQuadraturePoints(type);
  UInt nb_nodes_per_element = getFEM().getMesh().getNbNodesPerElement(type);

  std::stringstream name;
  name << id << ":structuralmechanics:imposed_linear_load";
  Vector<Real> traction_local(nb_element*nb_quad, nb_degree_of_freedom, name.str());

  Vector<Real>::const_iterator<types::RMatrix> T_it = rotation_matrix(type).begin(nb_degree_of_freedom * nb_nodes_per_element,
										  nb_degree_of_freedom * nb_nodes_per_element);

  Vector<Real>::const_iterator< types::RVector> Te_it = traction_global.begin(nb_degree_of_freedom);
  Vector<Real>::iterator< types::RVector> te_it = traction_local.begin(nb_degree_of_freedom);

  types::Matrix<Real> R(nb_degree_of_freedom, nb_degree_of_freedom);
  for (UInt e = 0; e < nb_element; ++e, ++T_it) {
    const types::RMatrix & T = *T_it;
    for (UInt i = 0; i < nb_degree_of_freedom; ++i)
      for (UInt j = 0; j < nb_degree_of_freedom; ++j)
	R(i, j) = T(i, j);

    for (UInt q = 0; q < nb_quad; ++q, ++Te_it, ++te_it) {
      const types::RVector & Te = *Te_it;
      types::RVector & te = *te_it;
      // turn the traction in the local referential
      te.mul<false>(R, Te);
    }
  }

  computeForcesByLocalTractionVector<type>(traction_local);

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
/**
 * @param myf pointer  to a function that fills a  vector/tensor with respect to
 * passed coordinates
 */
template<ElementType type>
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

  UInt offset = nb_degree_of_freedom;

  //prepare the loop over element types
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
  if(spatial_dimension == 3)
    getFEMClass<MyFEMType>().getShapeFunctions().interpolateOnControlPoints<type>(getFEM().getMesh().getNodes(),
										  quad_coords,
										  spatial_dimension,
										  _not_ghost,
										  NULL,
										  true,
										  0,
										  2,
										  2);
  lin_load.resize(nb_element*nb_quad);
  Real * imposed_val = lin_load.storage();

  /// sigma/load on each quadrature points
  Real * qcoord = quad_coords.storage();
  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quad; ++q) {
      myf(qcoord, imposed_val, NULL, 0);
      imposed_val += offset;
      qcoord += spatial_dimension;
    }
  }

  switch(function_type) {
  case _bft_traction_local:
    computeForcesByLocalTractionVector<type>(lin_load); break;
  case _bft_traction:
    computeForcesByGlobalTractionVector<type>(lin_load); break;
  default: break;
  }
}

/* -------------------------------------------------------------------------- */
