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
  registerParam("is_master_deformable", is_master_deformable, bool(false),
                _pat_parsable | _pat_readable, "Is master surface deformable");
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

  const auto slave_nodes =
      model.getContactDetector().getSurfaceSelector().getSlaveList();

  this->assembleInternalForces(slave_nodes);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleInternalForces(const Array<UInt> & slave_nodes) {
  AKANTU_DEBUG_IN();

  auto & contact_map = model.getContactMap();
    
  for (auto & slave : slave_nodes) {

    if (contact_map.find(slave) == contact_map.end())
      continue;

    auto & element    = contact_map[slave];
    const auto & conn = element.connectivity;
       
    Vector<Real> f_n(conn.size() * spatial_dimension);
    computeNormalForce(element, f_n);

    Vector<Real> f_t(conn.size() * spatial_dimension);
    computeTangentialForce(element, f_t);
 
    Vector<Real> f_c(conn.size() * spatial_dimension);
    f_c = f_n + f_t;

    assembleLocalToGlobalArray(slave, element, f_n, model.getNormalForce());
    assembleLocalToGlobalArray(slave, element, f_t, model.getTangentialForce());
    assembleLocalToGlobalArray(slave, element, f_c, model.getInternalForce());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleLocalToGlobalArray(const UInt & slave, const ContactElement & element,
					    Vector<Real> & local, Array<Real> & global) {

  auto get_connectivity = [&](auto & slave, auto & master) {
    Vector<UInt> master_conn =
    const_cast<const Mesh &>(model.getMesh()).getConnectivity(master);
    Vector<UInt> elem_conn(master_conn.size() + 1);

    elem_conn[0] = slave;
    for (UInt i = 1; i < elem_conn.size(); ++i) {
      elem_conn[i] = master_conn[i - 1];
    }

    return elem_conn;
  };

  auto connectivity = get_connectivity(slave, element.master);
  
  UInt nb_dofs  = global.getNbComponent();
  UInt nb_nodes = is_master_deformable ? connectivity.size() : 1;

  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());
  for (UInt i : arange(nb_nodes)) { 
    UInt n = connectivity[i];
    for (UInt j : arange(nb_dofs)) {
      UInt offset_node = n * nb_dofs + j;
      global[offset_node] += local[i * nb_dofs + j] * nodal_area[slave];
    }
  }
}
  
/* -------------------------------------------------------------------------- */
void Resolution::assembleStiffnessMatrix(GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  const auto slave_nodes =
      model.getContactDetector().getSurfaceSelector().getSlaveList();

  auto & stiffness =
      const_cast<SparseMatrix &>(model.getDOFManager().getMatrix("K"));

  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());

  auto & contact_map = model.getContactMap();

  for (auto & slave : slave_nodes) {

    if (contact_map.find(slave) == contact_map.end())
      continue;

    auto & element = contact_map[slave];

    const auto & conn = element.connectivity;

    Matrix<Real> kc(conn.size() * spatial_dimension,
                    conn.size() * spatial_dimension);

    Matrix<Real> m_alpha_beta(spatial_dimension - 1, spatial_dimension - 1);
    ResolutionUtils::computeMetricTensor(m_alpha_beta, element.tangents);

    // normal tangent moduli
    Vector<Real> n(conn.size() * spatial_dimension);
    ResolutionUtils::firstVariationNormalGap(element, n);
    
    Array<Real> t_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
    ResolutionUtils::computeTalpha(element, t_alpha);
    
    Array<Real> n_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
    ResolutionUtils::computeNalpha(element, n_alpha);
    
    Array<Real> d_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);      
    ResolutionUtils::firstVariationNaturalCoordinate(element, d_alpha);

    computeNormalModuli(kc, n_alpha, d_alpha, n, element);

    // frictional tangent moduli
    if (mu != 0) {
      /*Array<Real> t_alpha_beta(conn.size() * spatial_dimension,
                               (spatial_dimension - 1) * (spatial_dimension - 1));
      ResolutionUtils::computeTalphabeta(t_alpha_beta, element);
      
      Array<Real> p_alpha(conn.size() * spatial_dimension,
                          spatial_dimension - 1);
      Array<Real> n_alpha_beta(conn.size() * spatial_dimension,
                               (spatial_dimension - 1) * (spatial_dimension - 1));

      computeFrictionalTraction(m_alpha_beta, element);
      
      ResolutionUtils::computeNalphabeta(n_alpha_beta, element);
      ResolutionUtils::computePalpha(p_alpha, element);

      auto phi = computeNablaOfDisplacement(element);

      computeFrictionalModuli(kc, t_alpha_beta, n_alpha_beta, n_alpha, d_alpha,
      phi, n, element);*/
    }

    std::vector<UInt> equations;
    UInt nb_degree_of_freedom = model.getSpatialDimension();

    std::vector<Real> areas;
    
    UInt total_nodes = 1;
    UInt total_nb_degree_of_freedom = nb_degree_of_freedom;
    if (is_master_deformable) {
      total_nodes = conn.size();
      total_nb_degree_of_freedom *= total_nodes;
    }

    auto slave_node = conn[0];
    for (UInt i : arange(conn.size())) {
      UInt n = conn[i];
      for (UInt j : arange(nb_degree_of_freedom)) {
        equations.push_back(n * nb_degree_of_freedom + j);
        areas.push_back(nodal_area[slave_node]);
      }
    }

    for (UInt i : arange(total_nb_degree_of_freedom)) {
      UInt row = equations[i];
      for (UInt j : arange(total_nb_degree_of_freedom)) {
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

  // Matrix<Real> shape_second_derivatives(surface_dimension *
  // surface_dimension, 					nb_nodes_per_element);

  /*#define GET_SHAPE_SECOND_DERIVATIVES_NATURAL(type)			\
  ElementClass<type>::computeDN2DS2(element.projection, shape_second_derivatives)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_SECOND_DERIVATIVES_NATURAL);
  #undef GET_SHAPE_SECOND_DERIVATIVES_NATURAL*/

  Matrix<Real> nabla_u(surface_dimension * surface_dimension,
                       spatial_dimension);
  // nabla_u.mul<false, true>(shape_second_derivatives, values);

  return nabla_u;
}

} // namespace akantu
