/**
 * @file   structural_mechanics_model_inline_impl.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu May  5 19:48:07 2011
 *
 * @brief StructuralMechanicsModel implementation
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

/* -------------------------------------------------------------------------- */
template <ElementType type>
void StructuralMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  SparseMatrix & K = *stiffness_matrix;
  //  const Vector<Int> & equation_number = *equation_number;
  //  const Vector<Real> * shapes_derivatives;

  UInt nb_element                 = getFEM().getMesh().getNbElement(type);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(type);

  UInt tangent_size = getTangentStiffnessVoigtSize<type>();

  Vector<Real> * tangent_stiffness_matrix =
    new Vector<Real>(nb_element * nb_quadrature_points, tangent_size * tangent_size,
		     "tangent_stiffness_matrix");

  computeTangentStiffness<type>(*tangent_stiffness_matrix);

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

  Vector<Real>::iterator<types::RMatrix> B = b->begin(tangent_size, bt_d_b_size);
  Vector<Real>::iterator<types::RMatrix> D = tangent_stiffness_matrix->begin(tangent_size, tangent_size);
  Vector<Real>::iterator<types::RMatrix> Bt_D_B = bt_d_b->begin(bt_d_b_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Bt_D.mul<true, false>(*B, *D);
      Bt_D_B->mul<false, false>(Bt_D, *B);

      ++B;
      ++D;
      ++Bt_D_B;
    }
  }

  delete b;
  delete tangent_stiffness_matrix;

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
void StructuralMechanicsModel::computeTangentStiffness(__attribute__ ((unused)) Vector<Real> & tangent_stiffness_matrix) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


/* -------------------------------------------------------------------------- */
template<ElementType type>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix(__attribute__ ((unused)) Vector<Real> & b) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
void StructuralMechanicsModel::computeStressOnQuad() {
 AKANTU_DEBUG_IN();

 Vector<Real> & sigma  = stress(type, _not_ghost);

 sigma.clear();

  UInt nb_element                 = getFEM().getMesh().getNbElement(type);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = getFEM().getNbQuadraturePoints(type);
  const Vector<UInt> & connect    = getFEM().getMesh().getConnectivity(type);

  UInt tangent_size = getTangentStiffnessVoigtSize<type>();

  Vector<Real> * tangent_stiffness_matrix =
    new Vector<Real>(nb_element*nb_quadrature_points, tangent_size * tangent_size,
		     "tangent_stiffness_matrix");

  computeTangentStiffness<type>(*tangent_stiffness_matrix);

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
  Vector<Real>::iterator<types::RMatrix> D = tangent_stiffness_matrix->begin(tangent_size, tangent_size);
  Vector<Real>::iterator<types::RMatrix> D_B = d_b->begin(tangent_size, d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      D_B->mul<false, false>(*D, *B);

      ++B;
      ++D;
      ++D_B;
    }
  }

  delete b;
  delete tangent_stiffness_matrix;

  /// compute DBu
  D_B = d_b->begin(tangent_size, d_b_size);
  Vector<Real>::iterator< types::RVector> DBu = sigma.begin(tangent_size);
  types::RVector U (d_b_size);
  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      memcpy(U.storage()+n*nb_degree_of_freedom,
	     displacement_rotation->values+connect(e,n)*nb_degree_of_freedom,
	     nb_degree_of_freedom*sizeof(Real));
    }

    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      DBu ->mul<false>(*D_B,U);
      ++D_B;
      ++DBu;
    }
  }

  delete d_b;

  AKANTU_DEBUG_OUT();
}
