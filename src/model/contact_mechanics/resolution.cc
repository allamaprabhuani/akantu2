/**
 * @file   resolution.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Jan 7 2019
 * @date last modification: Mon Jan 7 2019
 *
 * @brief  Implementation of common part of the contact resolution class
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
#include "resolution.hh"
#include "contact_mechanics_model.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Resolution::Resolution(ContactMechanicsModel & model, const ID & id)
    : Memory(id, model.getMemoryID()),
      Parsable(ParserType::_contact_resolution, id), fem(model.getFEEngine()),
      name(""), model(model),
      spatial_dimension(model.getMesh().getSpatialDimension()) {

  AKANTU_DEBUG_IN();

  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Resolution::~Resolution() = default;

/* -------------------------------------------------------------------------- */
void Resolution::initialize() {
  registerParam("name", name, std::string(), _pat_parsable | _pat_readable);
  registerParam("mu", mu, Real(0.), _pat_parsable | _pat_modifiable,
                "Friciton Coefficient");
}

/* -------------------------------------------------------------------------- */
void Resolution::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  std::string type = getID().substr(getID().find_last_of(':') + 1);

  stream << space << "Contact Resolution " << type << " [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleInternalForces(GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  auto & internal_force = const_cast<Array<Real> &>(model.getInternalForce());

  const auto local_nodes = model.getMesh().getElementGroup(name).getNodes();

  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());

  auto & contact_map = model.getContactMap();

  for (auto & slave : local_nodes) {

    if (contact_map.find(slave) == contact_map.end())
      continue;

    auto & element = contact_map[slave];

    const auto & conn = element.connectivity;
    const auto & type = element.master.type;

    auto nb_nodes_master = Mesh::getNbNodesPerElement(type);

    Vector<Real> shapes(nb_nodes_master);
    Matrix<Real> dnds(spatial_dimension - 1, nb_nodes_master);

#define GET_SHAPES_NATURAL(type)                                               \
    ElementClass<type>::computeShapes(element.projection, shapes)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL

#define GET_SHAPE_DERIVATIVES_NATURAL(type)                                    \
    ElementClass<type>::computeDNDS(element.projection, dnds)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL

    Vector<Real> fc(conn.size() * spatial_dimension);
   
    Vector<Real> n(conn.size() * spatial_dimension);
    computeN(n, shapes, element.normal);
    
    computeNormalForce(fc, n, element.gap);

    if(mu != 0) {
      Matrix<Real> m_alpha_beta(spatial_dimension - 1, spatial_dimension - 1);
      computeMetricTensor(element.tangents, m_alpha_beta);

      Array<Real> t_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
      Array<Real> n_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
      Array<Real> d_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);

      computeTalpha(t_alpha, shapes,
		    element.tangents);
      computeNalpha(n_alpha, dnds,
		    element.normal);
      computeDalpha(d_alpha, n_alpha, t_alpha,
		    m_alpha_beta, element.gap);
    
      auto traction = computeFrictionalTraction(m_alpha_beta,
						element.projection,
						element.gap);
      computeFrictionForce(fc, d_alpha, traction);
    }
    
    UInt nb_degree_of_freedom = internal_force.getNbComponent();
    for (UInt i = 0; i < conn.size(); ++i) {

      UInt n = conn[i];
      for (UInt j = 0; j < nb_degree_of_freedom; ++j) {
        UInt offset_node = n * nb_degree_of_freedom + j;
        internal_force[offset_node] += fc[i * nb_degree_of_freedom + j];
        internal_force[offset_node] *= nodal_area[n];
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleStiffnessMatrix(GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  auto & stiffness =
      const_cast<SparseMatrix &>(model.getDOFManager().getMatrix("K"));

  const auto local_nodes =
    model.getMesh().getElementGroup(name).getNodes();

  auto & nodal_area =
    const_cast<Array<Real> &>(model.getNodalArea());

  auto & contact_map = model.getContactMap();

  for (auto & slave : local_nodes) {

    if (contact_map.find(slave) == contact_map.end())
      continue;
    
    auto & element = contact_map[slave];
                
    const auto & conn = element.connectivity;
    const auto & type = element.master.type;
    
    UInt nb_nodes_master = Mesh::getNbNodesPerElement(type);

    Vector<Real> shapes(nb_nodes_master);
    Matrix<Real> shapes_derivatives(spatial_dimension - 1, nb_nodes_master);
    Matrix<Real> shapes_second_derivatives((spatial_dimension-1)*(spatial_dimension-1) ,
					   nb_nodes_master);
    
#define GET_SHAPES_NATURAL(type)					\
    ElementClass<type>::computeShapes(element.projection, shapes)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL

#define GET_SHAPE_DERIVATIVES_NATURAL(type)				\
    ElementClass<type>::computeDNDS(element.projection, shapes_derivatives)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL

#define GET_SHAPE_SECOND_DERIVATIVES_NATURAL(type)			\
    ElementClass<type>::computeDN2DS2(element.projection, shapes_second_derivatives)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_SECOND_DERIVATIVES_NATURAL);
#undef GET_SHAPE_SECOND_DERIVATIVES_NATURAL

    Matrix<Real> kc(conn.size() * spatial_dimension,
		    conn.size() * spatial_dimension);
    
    Matrix<Real> m_alpha_beta(spatial_dimension - 1, spatial_dimension - 1);
    computeMetricTensor(element.tangents, m_alpha_beta);

    // normal tangent moduli
    Vector<Real> n(conn.size() * spatial_dimension);
    Array<Real> t_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
    Array<Real> n_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
    Array<Real> d_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);

    computeN(      n,        shapes,             element.normal);
    computeTalpha( t_alpha,  shapes,             element.tangents);
    computeNalpha( n_alpha,  shapes_derivatives, element.normal);
    computeDalpha( d_alpha,  n_alpha,  t_alpha,  m_alpha_beta, element.gap);

    computeNormalModuli(kc, n, n_alpha, d_alpha, m_alpha_beta, element.gap);

    // frictional tangent moduli
    if(mu != 0) {
      Array<Real> t_alpha_beta(conn.size() * spatial_dimension,
			       (spatial_dimension - 1) * (spatial_dimension -1));
      Array<Real> p_alpha(conn.size() * spatial_dimension,
			  spatial_dimension - 1);
      Array<Real> n_alpha_beta(conn.size() * spatial_dimension,
			       (spatial_dimension - 1) * (spatial_dimension -1));

      auto traction = computeFrictionalTraction(m_alpha_beta, element.projection,
						element.gap);
      computeTalphabeta(t_alpha_beta, shapes_derivatives, element.tangents);
      computeNalphabeta(n_alpha_beta, shapes_second_derivatives, element.normal);
      computePalpha(p_alpha, shapes_derivatives, traction);
      //computeGalpha();
      //computePbaralpha();
      //computeTbaralphabeta();
      
      computeFrictionalModuli(kc, t_alpha_beta, n_alpha_beta, element.tangents,
			      shapes_second_derivatives, n, n_alpha, d_alpha,
			      element.gap);
    }
        
    std::vector<UInt> equations;
    UInt nb_degree_of_freedom = model.getSpatialDimension();

    std::vector<Real> areas;
    for (UInt i : arange(connectivity.size())) {
      UInt n = connectivity[i];
      for (UInt j : arange(nb_degree_of_freedom)) {
        equations.push_back(n * nb_degree_of_freedom + j);
        areas.push_back(nodal_area[n]);
      }
    }

    for (UInt i : arange(kc.rows())) {
      UInt row = equations[i];
      for (UInt j : arange(kc.cols())) {
        UInt col = equations[j];
        kc(i, j) *= areas[i];
        stiffness.add(row, col, kc(i, j));
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */ 
void Resolution::computeMetricTensor(Matrix<Real> & tangents, Matrix<Real> & m_alpha_beta) {

  m_alpha_beta.mul<false, true>(tangents, tangents);
  m_alpha_beta = m_alpha_beta.inverse();
}

/* -------------------------------------------------------------------------- */
void Resolution::computeN(Vector<Real> & n, Vector<Real> & shapes,
                          Vector<Real> & normal) {

  UInt dim = normal.size();
  for (UInt i = 0; i < dim; ++i) {
    n[i] = normal[i];
    for (UInt j = 0; j < shapes.size(); ++j) {
      n[(1 + j) * dim + i] = -normal[i] * shapes[j];
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeTalpha(Array<Real> & t_alpha, Vector<Real> & shapes,
                               Matrix<Real> & tangents) {
  t_alpha.clear();
  for (auto && values :
       zip(tangents.transpose(), make_view(t_alpha, t_alpha.size()))) {

    auto & tangent = std::get<0>(values);
    auto & t_s = std::get<1>(values);
    for (UInt i : arange(spatial_dimension)) {
      t_s[i] = -tangent(i);
      for (UInt j : arange(shapes.size())) {
        t_s[(1 + j) * spatial_dimension + i] = -shapes[j] * tangent(i);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeNalpha(Array<Real> & n_alpha,
                               Matrix<Real> & shapes_derivatives,
                               Vector<Real> & normal) {
  n_alpha.clear();
  for (auto && values : zip(shapes_derivatives.transpose(),
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
void Resolution::computeDalpha(Array<Real> & d_alpha, Array<Real> & n_alpha,
			       Array<Real> & t_alpha, Matrix<Real> & m_alpha_beta,
			       Real & gap) {

  d_alpha.clear();
  for (auto && entry : zip(m_alpha_beta.transpose(),
			   make_view(d_alpha, d_alpha.size()))) {

    auto & a_s = std::get<0>(entry);
    auto & d_s = std::get<1>(entry);
    for (auto && values :
         zip(arange(t_alpha.size()),
	     make_view(t_alpha, t_alpha.size()),
             make_view(n_alpha, n_alpha.size()))) {
      auto & index = std::get<0>(values);
      auto & t_s = std::get<1>(values);
      auto & n_s = std::get<2>(values);

      d_s += (t_s + gap * n_s);
      d_s *= a_s(index);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeTalphabeta(Array<Real> & t_alpha_beta,
				   Matrix<Real> & shapes_derivatives,
				   Matrix<Real> & tangents) {
  t_alpha_beta.clear();

  auto t_alpha_size = t_alpha_beta.size() * (spatial_dimension - 1);
  for(auto && entry : zip(tangents.transpose(),
			  make_view(t_alpha_beta, t_alpha_size))) {
    auto & tangent_s = std::get<0>(entry);
    auto & t_alpha   = std::get<1>(entry);
    for(auto && values : zip(shapes_derivatives.transpose(),
    			     make_view(t_alpha, t_alpha_beta.size()))) {
      auto & dnds      = std::get<0>(values);
      auto & t_alpha_s = std::get<1>(values);

      for (UInt i : arange(spatial_dimension)) {
	t_alpha_s[i] = 0;
	for (UInt j : arange(dnds.size())) {
	  t_alpha_s[(1 + j) * spatial_dimension + i] = -dnds(j) * tangent_s[i];
	}
      } 
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeNalphabeta(Array<Real> & n_alpha_beta,
				   Matrix<Real> & shapes_second_derivatives,
				   Vector<Real> & normal) {
  n_alpha_beta.clear();

  for(auto && entry : zip(shapes_second_derivatives.transpose(),
			  make_view(n_alpha_beta, n_alpha_beta.size()))) {
    auto & dn2ds2    = std::get<0>(entry);
    auto & n_alpha_s = std::get<1>(entry);

    for (UInt i : arange(spatial_dimension)) {
      n_alpha_s[i] = 0;
      for (UInt j : arange(dn2ds2.size())) {
        n_alpha_s[(1 + j) * spatial_dimension + i] = -dn2ds2(j) * normal[i];
      }
    }
       
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computePalpha(Array<Real> & p_alpha,
			       Matrix<Real> & shapes_derivatives,
			       Vector<Real> & traction) {
  p_alpha.clear();

  auto normalized_traction = traction/traction.norm();
  
  for(auto && entry :
	zip(shapes_derivatives.transpose(),
	    make_view(p_alpha, p_alpha.size()))) {
    auto & dnds = std::get<0>(entry);
    auto & p_s  = std::get<1>(entry);

    for(UInt i : arnage(spatial_dimension)) {
      p_s[i] = 0;
      for(UInt j : arange(dnds.size())){
	p_s[(1 + j) * spatial_dimension + i] = -dnds(j) * normalized_traction[i];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeGalpha(Array<Real> & /*t_alpha_beta*/,
			       Array<Real> & /*d_alpha*/,
			       Vector<Real> & /*tangential_gap*/) {

  Array<Real> g_alpha(d_alpha.size(), spatial_dimension - 1);
  auto t_alpha_size = t_alpha_beta.size() * (spatial_dimension - 1);
  
  for(auto && value :
	zip(make_view(g_alpha, g_alpha.size()),
	    make_view(t_alpha_beta, t_alpha_size))){
    auto & g_s = std::get<0>(value);
    
  }  
}

/* -------------------------------------------------------------------------- */
void Resolution::computeTbaralphabeta()  {

}

/* -------------------------------------------------------------------------- */
void Resolution::computePbaralpha() {

}
  
/* -------------------------------------------------------------------------- */
void Resolution::computeCoordinates(const Element & el, Matrix<Real> & coords) {
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
  Vector<UInt> connect = model.getMesh()
                             .getConnectivity(el.type, _not_ghost)
                             .begin(nb_nodes_per_element)[el.element];

  // todo change this to current position
  auto & positions = model.getMesh().getNodes();
  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt node = connect[n];
    for (UInt s : arange(spatial_dimension)) {
      coords(n, s) = positions(node, s);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeSecondDerivative(Matrix<Real> & dn2ds2, const Element & el) {

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
  Vector<UInt> connect = model.getMesh().getConnectivity(el.type, _not_ghost).begin(nb_nodes_per_element)[el.element];
  
  Matrix<Real> phi(dn2ds2.rows(), spatial_dimension);
  Matrix<Real> nodal_values(nb_nodes_per_element, spatial_dimension);

  computeCoordinates(el, nodal_values);
  phi.mul<false, false>(dn2ds2, nodal_values);
}
  
} // namespace akantu
