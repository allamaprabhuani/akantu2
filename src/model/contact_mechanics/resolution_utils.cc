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
void ResolutionUtils::firstVariationNormalGap(const ContactElement & element,
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
}

/* -------------------------------------------------------------------------- */
void ResolutionUtils::secondVariationNormalGap(const ContactElement & element,
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
}
  
/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeTalpha(const ContactElement & element,
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
}

/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeNalpha(const ContactElement & element, const Vector<Real> & projection,
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
}

/* -------------------------------------------------------------------------- */
void ResolutionUtils::firstVariationNaturalCoordinate(const ContactElement & element,
						      const Matrix<Real> & covariant_basis,
						      const Vector<Real> & projection,
						      const Vector<Real> & normal, const Real & gap,
						      Array<Real> & delta_xi) {

  delta_xi.clear();
  
  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;
  
  Matrix<Real> m_alpha_beta(surface_dimension, surface_dimension);
  ResolutionUtils::computeMetricTensor(m_alpha_beta, covariant_basis);
  m_alpha_beta = m_alpha_beta.inverse();

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
      auto & n_beta = std::get<2>(values);

      d_alpha += (t_beta + gap * n_beta) * m_alpha_beta(alpha, beta);
    }
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeTalphabeta(Array<Real> & t_alpha_beta,
					ContactElement & element) {
  t_alpha_beta.clear();
  
  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
    
  Matrix<Real> shape_derivatives(surface_dimension,
				 Mesh::getNbNodesPerElement(type));

#define GET_SHAPE_DERIVATIVES_NATURAL(type)				\
  ElementClass<type>::computeDNDS(element.projection, shape_derivatives)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL

  //auto t_alpha_size = t_alpha_beta.size() * surface_dimension;

  //auto & tangents = element.tangents;
  /*for(auto && entry :
	zip(tangents.transpose(),
	    make_view(t_alpha_beta, t_alpha_size))) {

    auto & tangent_s = std::get<0>(entry);
    auto & t_alpha   = std::get<1>(entry);

    for(auto && values :
	  zip(shape_derivatives.transpose(),
	      make_view(t_alpha, t_alpha_beta.size()))) {

      auto & dnds      = std::get<0>(values);
      auto & t_alpha_s = std::get<1>(values);

      for (UInt i : arange(spatial_dimension)) {
	t_alpha_s[i] = 0;
	for (UInt j : arange(dnds.size())) {
	  t_alpha_s[(1 + j) * spatial_dimension + i] = -dnds(j) * tangent_s(i);
	}
      } 
    }
  }*/
}

/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeNalphabeta(Array<Real> & n_alpha_beta,
					ContactElement & element) {
  n_alpha_beta.clear();

  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;

  Matrix<Real> shape_second_derivatives(surface_dimension * surface_dimension,
					Mesh::getNbNodesPerElement(type));
  
#define GET_SHAPE_SECOND_DERIVATIVES_NATURAL(type)			\
    ElementClass<type>::computeDN2DS2(element.projection, shape_second_derivatives)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_SECOND_DERIVATIVES_NATURAL);
#undef GET_SHAPE_SECOND_DERIVATIVES_NATURAL

  for(auto && entry :
	zip(shape_second_derivatives.transpose(),
	    make_view(n_alpha_beta, n_alpha_beta.size()))) {

    auto & dn2ds2    = std::get<0>(entry);
    auto & n_alpha_s = std::get<1>(entry);

    for (UInt i : arange(spatial_dimension)) {
      n_alpha_s[i] = 0;
      for (UInt j : arange(dn2ds2.size())) {
        n_alpha_s[(1 + j) * spatial_dimension + i] = -dn2ds2(j) * element.normal[i];
      }
    }
       
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionUtils::computePalpha(Array<Real> & p_alpha,
				    ContactElement & /*element*/) {
  p_alpha.clear();

  /*const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;
  
  Matrix<Real> shape_derivatives(surface_dimension,
				 Mesh::getNbNodesPerElement(type));

#define GET_SHAPE_DERIVATIVES_NATURAL(type)				\
  ElementClass<type>::computeDNDS(element.projection, shape_derivatives)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL

  auto normalized_traction = element.traction/element.traction.norm();
  
  for(auto && entry :
	zip(shape_derivatives.transpose(),
	    make_view(p_alpha, p_alpha.size()))) {
    auto & dnds = std::get<0>(entry);
    auto & p_s  = std::get<1>(entry);

    for(UInt i : arange(spatial_dimension)) {
      p_s[i] = 0;
      for(UInt j : arange(dnds.size())){
	p_s[(1 + j) * spatial_dimension + i] = -dnds(j) * normalized_traction[i];
      }
    }
  }*/
}

/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeGalpha(Array<Real> & g_alpha, Array<Real> & t_alpha_beta,
				    Array<Real> & d_alpha, Matrix<Real> & phi,
				    ContactElement & element) {

  g_alpha.clear();

  const auto & type = element.master.type;
  auto surface_dimension = Mesh::getSpatialDimension(type);

  auto & tangents = element.tangents;
  auto tangential_gap = element.projection - element.previous_projection;

  for(auto alpha : arange(surface_dimension)) {
    auto & g_a = g_alpha[alpha];
    auto & tangent_alpha = tangents[alpha];
    
    for(auto beta : arange(surface_dimension)) {
      auto & t_a_b = t_alpha_beta(alpha, beta); 
      auto & t_b_a = t_alpha_beta(beta, alpha);
      auto & gt_beta = tangential_gap[beta];
      auto & tangents_beta = tangents[beta];
      
      for(auto gamma : arange(surface_dimension)) {
	auto & d_gamma = d_alpha(gamma);
	auto tmp = phi(beta, gamma) * tangent_alpha + phi(alpha,gamma) * tangents_beta;

	g_a += (-t_a_b - t_b_a + tmp * d_gamma)*gt_beta;
      }
    }

  }

}

/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeMetricTensor(Matrix<Real> & m_alpha_beta, const Matrix<Real> & tangents) {

  m_alpha_beta.mul<false, true>(tangents, tangents);
}


/* -------------------------------------------------------------------------- */
void ResolutionUtils::computeSecondMetricTensor(const ContactElement & element,
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
}
  
/* -------------------------------------------------------------------------- */
void ResolutionUtils::assembleToInternalForce(Vector<Real> & local_array,
                                              Array<Real> & global_array,
                                              Array<Real> & nodal_area,
                                              ContactElement & element, bool is_master_deformable) {
  const auto & conn = element.connectivity;
  UInt nb_dofs = global_array.getNbComponent();

  auto slave_node = conn[0];
 
  UInt total_nodes = 1;
  if (is_master_deformable) {
    total_nodes = conn.size();
  }
  
  for (UInt i : arange(total_nodes)) { // 1 to only consider slave node 
    UInt n = conn[i];
    for (UInt j : arange(nb_dofs)) {
      UInt offset_node = n * nb_dofs + j;
      global_array[offset_node] += local_array[i * nb_dofs + j] * nodal_area[slave_node];
    }
  }
}

/* -------------------------------------------------------------------------- */
  void ResolutionUtils::assembleToStiffnessMatrix(Matrix<Real> & /*local_matrix*/, Matrix<Real> & /*global_matrix*/,
						  Array<Real> & /*nodal_area*/, ContactElement & /*element*/) {
  
}
  
  
} //akantu

