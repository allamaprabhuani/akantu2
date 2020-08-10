/**
 * @file   resolution_utils.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Mmay 20 2019
 * @date last modification: Mon May 20 2019
 *
 * @brief  Implementation of various utilities neede for resolution class 
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "resolution_utils.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeShapeFunctionMatric(const ContactElement & element,
						 const Vector<Real> & projection,
						 Matrix<Real> & shape_matric) {
  
  shape_matric.clear();

  const ElementType & type = element.master.type;
  
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;
  UInt nb_nodes_per_contact = element.getNbNodes();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  AKANTU_DEBUG_ASSERT(
	spatial_dimension == shape_matric.rows()
	&& spatial_dimension* nb_nodes_per_contact == shape_matric.cols(),
        "Shape Matric dimensions are not correct");

  Vector<Real> shapes(nb_nodes_per_element);  

#define GET_SHAPE_NATURAL(type)				\
  ElementClass<type>::computeShapes(projection, shapes)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_NATURAL);
#undef GET_SHAPE_NATURAL

  for (auto i : arange(nb_nodes_per_contact)) {
    for (auto j : arange(spatial_dimension)) {
      if (i == 0) {
	shape_matric(j, i*spatial_dimension + j) = 1;
	continue;
      }
      shape_matric(j, i*spatial_dimension + j) = -shapes[i-1];
    }
  }

}

/* -------------------------------------------------------------------------- */
/*void ResolutionUtils::firstVariationNormalGap(const ContactElement & element,
					      const Vector<Real> & projection,
					      const Vector<Real> & normal,
					      Vector<Real> & delta_g) {

  delta_g.clear();
  
  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;
  
  Vector<Real> shapes(Mesh::getNbNodesPerElement(type));
  
#define GET_SHAPES_NATURAL(type)                                               \
    ElementClass<type>::computeShapes(projection, shapes)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL

  for (UInt i : arange(spatial_dimension)) {
    delta_g[i] = normal[i];
    for (UInt j : arange(shapes.size())) {
      delta_g[(1 + j) * spatial_dimension + i] = -shapes[j] * normal[i];
    }
  }
  }*/

/* -------------------------------------------------------------------------- */
/*void ResolutionUtils::secondVariationNormalGap(const ContactElement & element,
					       const Matrix<Real> & covariant_basis,
					       const Matrix<Real> & curvature,
 					       const Vector<Real> & projection,
					       const Vector<Real> & normal, Real & gap,
					       Matrix<Real> & ddelta_g) {

  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;

  UInt nb_nodes = element.getNbNodes();
  
  Array<Real> dnds_n(nb_nodes * spatial_dimension, surface_dimension);
  ResolutionUtils::computeNalpha(element, projection, normal, dnds_n);
  
  Array<Real> delta_xi(nb_nodes * spatial_dimension, surface_dimension);
  ResolutionUtils::firstVariationNaturalCoordinate(element, covariant_basis,
						   projection, normal, gap, delta_xi);
  
  Matrix<Real> a_alpha_beta(surface_dimension, surface_dimension);
  ResolutionUtils::computeMetricTensor(a_alpha_beta, covariant_basis);
  a_alpha_beta = a_alpha_beta.inverse();

  Matrix<Real> h_alpha_beta(surface_dimension, surface_dimension);
  ResolutionUtils::computeSecondMetricTensor(element, curvature,
					     normal, h_alpha_beta);

  for (auto && values : zip(arange(surface_dimension),
			    make_view(dnds_n, dnds_n.size()),
			    make_view(delta_xi, delta_xi.size()))) {

    auto & alpha = std::get<0>(values);
    auto & dnds_n_alpha = std::get<1>(values);
    auto & delta_xi_alpha = std::get<2>(values);

    // term 1 from Numerical methods in contact mechanics : Vlad
    // Yastrebov eq 2.48 
    Matrix<Real> mat_n(dnds_n_alpha.storage(), dnds_n_alpha.size(), 1);
    Matrix<Real> mat_xi(delta_xi_alpha.storage(), delta_xi_alpha.size(), 1);

    Matrix<Real> tmp1(dnds_n_alpha.size(), dnds_n_alpha.size());
    tmp1.mul<false, true>(mat_n, mat_xi, -1);

    Matrix<Real> tmp2(dnds_n_alpha.size(), dnds_n_alpha.size());
    tmp2.mul<false, true>(mat_xi, mat_n, -1);

    Matrix<Real> term1(dnds_n_alpha.size(), dnds_n_alpha.size());
    term1 = tmp1 + tmp2;

    // computing term 2 & term 3 from Numerical methods in contact
    // mechanics : Vlad Yastrebov eq 2.48
    Matrix<Real> term2(delta_xi_alpha.size(), delta_xi_alpha.size());
    Matrix<Real> term3(dnds_n_alpha.size(), dnds_n_alpha.size());

    for (auto && values2 : zip(arange(surface_dimension),
			       make_view(dnds_n, dnds_n.size()),
			       make_view(delta_xi, delta_xi.size()))) {
      auto & beta = std::get<0>(values2);
      auto & dnds_n_beta = std::get<1>(values2);
      auto & delta_xi_beta = std::get<2>(values2);

      // term 2
      Matrix<Real> mat_xi_beta(delta_xi_beta.storage(), delta_xi.size(), 1);
      Matrix<Real> tmp3(delta_xi_beta.size(), delta_xi_beta.size());

      Real pre_factor = h_alpha_beta(alpha, beta);
      for (auto k : arange(surface_dimension)) {
	for (auto m : arange(surface_dimension)) {
	  pre_factor -= gap * h_alpha_beta(alpha, k) * a_alpha_beta(k, m) * h_alpha_beta(m, beta);
	}
      }

      pre_factor *= -1.;
      tmp3.mul<false, true>(mat_xi, mat_xi_beta, pre_factor);

      // term 3
      Matrix<Real> mat_n_beta(dnds_n_beta.storage(), dnds_n_beta.size(), 1);

      Real factor = gap * a_alpha_beta(alpha, beta);
      Matrix<Real> tmp4(dnds_n_alpha.size(), dnds_n_alpha.size());
      tmp4.mul<false, true>(mat_n, mat_n_beta, factor);

      term3 += tmp4;
    }

    ddelta_g += term1 + term2 + term3;
  }
  }*/
  
/* -------------------------------------------------------------------------- */
/*void ResolutionUtils::computeTalpha(const ContactElement & element,
				    const Matrix<Real> & covariant_basis,
				    const Vector<Real> & projection, Array<Real> & t_alpha) {

  t_alpha.clear();

  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;
  
  Vector<Real> shapes(Mesh::getNbNodesPerElement(type));

#define GET_SHAPES_NATURAL(type)                                               \
  ElementClass<type>::computeShapes(projection, shapes)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL

  for (auto && values :
	 zip(covariant_basis.transpose(),
	     make_view(t_alpha, t_alpha.size()))) {

    auto & tangent_beta = std::get<0>(values);
    auto & t_beta       = std::get<1>(values);
    
    for (UInt i : arange(spatial_dimension)) {
      t_beta[i] = tangent_beta(i);
      for (UInt j : arange(shapes.size())) {
	t_beta[(1 + j) * spatial_dimension + i] = -shapes[j] * tangent_beta(i);
      }
    }
  }
  }*/

/* -------------------------------------------------------------------------- */
/*void ResolutionUtils::computeNalpha(const ContactElement & element, const Vector<Real> & projection,
				    const Vector<Real> & normal, Array<Real> & n_alpha) {

  n_alpha.clear();

  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;
  
  Matrix<Real> shape_derivatives(surface_dimension,
				 Mesh::getNbNodesPerElement(type));

#define GET_SHAPE_DERIVATIVES_NATURAL(type)				\
  ElementClass<type>::computeDNDS(projection, shape_derivatives)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL
  
  for (auto && values :
	 zip(shape_derivatives.transpose(),
	     make_view(n_alpha, n_alpha.size()))) {

    auto & dnds = std::get<0>(values);
    auto & n_s = std::get<1>(values);

    for (UInt i : arange(spatial_dimension)) {
      n_s[i] = 0;
      for (UInt j : arange(dnds.size())) {
        n_s[(1 + j) * spatial_dimension + i] = -dnds(j) * normal[i];
      }
    }
  }
  }*/

/* -------------------------------------------------------------------------- */
/*void ResolutionUtils::firstVariationNaturalCoordinate(const ContactElement & element,
						      const Matrix<Real> & covariant_basis,
						      const Vector<Real> & projection,
						      const Vector<Real> & normal, const Real & gap,
						      Array<Real> & delta_xi) {

  delta_xi.clear();
  
  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;
  
  auto inv_A = GeometryUtils::contravariantMetricTensor(covariant_basis);

  auto nb_nodes = element.getNbNodes();
  
  Array<Real> t_alpha(nb_nodes * spatial_dimension, surface_dimension);
  Array<Real> n_alpha(nb_nodes * spatial_dimension, surface_dimension);

  ResolutionUtils::computeTalpha(element, covariant_basis, projection, t_alpha);
  ResolutionUtils::computeNalpha(element, projection, normal, n_alpha);
  
  for (auto && entry :
	 zip(arange(surface_dimension),
	     make_view(delta_xi, delta_xi.size()))) {

    auto & alpha   = std::get<0>(entry);
    auto & d_alpha = std::get<1>(entry);

    for (auto && values :
         zip(arange(surface_dimension),
	     make_view(t_alpha, t_alpha.size()),
             make_view(n_alpha, n_alpha.size()))) {
      auto & beta   = std::get<0>(values);
      auto & t_beta = std::get<1>(values);
      //auto & n_beta = std::get<2>(values);

      //d_alpha += (t_beta + gap * n_beta) * m_alpha_beta(alpha,
      //beta);
      d_alpha += t_beta * inv_A(alpha, beta);
    }
  }
}*/

/* -------------------------------------------------------------------------- */
/*void ResolutionUtils::computeMetricTensor(Matrix<Real> & m_alpha_beta,
					  const Matrix<Real> & tangents) {

  m_alpha_beta.mul<false, true>(tangents, tangents);
  }*/


/* -------------------------------------------------------------------------- */
/*void ResolutionUtils::computeSecondMetricTensor(const ContactElement & element,
						const Matrix<Real> & curvature,
						const Vector<Real> & normal,
						Matrix<Real> & metric) {

  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  
  auto i = 0;
  for (auto alpha : arange(surface_dimension) ) {
    for (auto beta : arange(surface_dimension)) {
      Vector<Real> temp(curvature(i));
      metric(alpha, beta) = normal.dot(temp);
      i++;
    }
  }    
  }*/
  
  
  
} //akantu

