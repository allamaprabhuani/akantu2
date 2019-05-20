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

    Vector<Real> contact_force(conn.size() * spatial_dimension);
   
    Vector<Real> n(conn.size() * spatial_dimension);
    ResolutionUtils::computeN(n, element);
    
    computeNormalForce(contact_force, n, element);

    if(mu != 0) {

      Array<Real> t_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
      Array<Real> n_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
      Array<Real> d_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);

      ResolutionUtils::computeTalpha(t_alpha, element);
      ResolutionUtils::computeNalpha(n_alpha, element);
      ResolutionUtils::computeDalpha(d_alpha, n_alpha, t_alpha, element);
   
      computeFrictionalForce(contact_force, d_alpha, element);
    }
    
    ResolutionUtils::assembleToInternalForce(contact_force, internal_force,
					     nodal_area, element);
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

    Matrix<Real> kc(conn.size() * spatial_dimension,
		    conn.size() * spatial_dimension);
    
    Matrix<Real> m_alpha_beta(spatial_dimension - 1, spatial_dimension - 1);
    ResolutionUtils::computeMetricTensor(element.tangents, m_alpha_beta);

    // normal tangent moduli
    Vector<Real> n(conn.size() * spatial_dimension);

    Array<Real> t_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
    Array<Real> n_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
    Array<Real> d_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);

    ResolutionUtils::computeN(      n,        element);
    ResolutionUtils::computeTalpha( t_alpha,  element);
    ResolutionUtils::computeNalpha( n_alpha,  element);
    ResolutionUtils::computeDalpha( d_alpha,  n_alpha,  t_alpha, element);

    computeNormalModuli(kc, n_alpha, d_alpha, n, element);

    // frictional tangent moduli
    if(mu != 0) {
      Array<Real> t_alpha_beta(conn.size() * spatial_dimension,
			       (spatial_dimension - 1) * (spatial_dimension -1));
      Array<Real> p_alpha(conn.size() * spatial_dimension,
			  spatial_dimension - 1);
      Array<Real> n_alpha_beta(conn.size() * spatial_dimension,
			       (spatial_dimension - 1) * (spatial_dimension -1));

      computeFrictionalTraction(m_alpha_beta, element);

      ResolutionUtils::computeTalphabeta(t_alpha_beta, element);
      ResolutionUtils::computeNalphabeta(n_alpha_beta, element);
      ResolutionUtils::computePalpha(p_alpha, element);

      auto phi = computeNablaOfDisplacement(element);
            
      computeFrictionalModuli(kc, t_alpha_beta, n_alpha_beta,
			      n_alpha, d_alpha, phi, n, element);
    }
        
    std::vector<UInt> equations;
    UInt nb_degree_of_freedom = model.getSpatialDimension();

    std::vector<Real> areas;
    for (UInt i : arange(conn.size())) {
      UInt n = conn[i];
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
Matrix<Real> Resolution::computeNablaOfDisplacement(ContactElement & element) {

  const auto & type = element.master.type;
  const auto & conn = element.connectivity;
  
  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;

  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  Matrix<Real> values(spatial_dimension, nb_nodes_per_element);
  
  auto & displacement = model.getDisplacement();
  for (UInt n : arange(nb_nodes_per_element)) {
    UInt node = conn[n];
    for (UInt s : arange(spatial_dimension)) {
      values(s, n) = displacement(node, s);
    }
  }

  Matrix<Real> shape_second_derivatives(surface_dimension * surface_dimension,
					nb_nodes_per_element);
  
  //#define GET_SHAPE_SECOND_DERIVATIVES_NATURAL(type)			\
  //ElementClass<type>::computeDN2DS2(element.projection, shape_second_derivatives)
  // AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_SECOND_DERIVATIVES_NATURAL);
  //#undef GET_SHAPE_SECOND_DERIVATIVES_NATURAL

  Matrix<Real> nabla_u(surface_dimension * surface_dimension, spatial_dimension);
  //nabla_u.mul<false, true>(shape_second_derivatives, values);

  return nabla_u;
}
  
} // namespace akantu
