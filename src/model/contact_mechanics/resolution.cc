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
      name(""), model(model) {

  AKANTU_DEBUG_IN();

  spatial_dimension = model.getMesh().getSpatialDimension();
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

  this->assembleInternalForces();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  for (auto & element : model.getContactElements()) {
     
    auto nb_nodes  = element.getNbNodes();
           
    Vector<Real> f_n(nb_nodes * spatial_dimension);
    computeNormalForce(element, f_n);

    Vector<Real> f_t(nb_nodes * spatial_dimension);
    computeTangentialForce(element, f_t);
 
    Vector<Real> f_c(nb_nodes * spatial_dimension);
    f_c = f_n + f_t;

    assembleLocalToGlobalArray(element, f_n, model.getNormalForce());
    assembleLocalToGlobalArray(element, f_t, model.getTangentialForce());
    assembleLocalToGlobalArray(element, f_c, model.getInternalForce());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleLocalToGlobalArray(const ContactElement & element,
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

  auto connectivity = get_connectivity(element.slave, element.master);
  
  UInt nb_dofs  = global.getNbComponent();
  UInt nb_nodes = is_master_deformable ? connectivity.size() : 1;

  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());
  for (UInt i : arange(nb_nodes)) { 
    UInt n = connectivity[i];
    for (UInt j : arange(nb_dofs)) {
      UInt offset_node = n * nb_dofs + j;
      global[offset_node] += local[i * nb_dofs + j] * nodal_area[element.slave];
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleLocalToGlobalMatrix(const ContactElement & element,
					     const Matrix<Real> & local, SparseMatrix & global) {

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

  auto connectivity = get_connectivity(element.slave, element.master);

  auto nb_dofs  = spatial_dimension;
  UInt nb_nodes = is_master_deformable ? connectivity.size() : 1;
  UInt total_nb_dofs = nb_dofs * nb_nodes;
  
  std::vector<UInt> equations;     
  for (UInt i : arange(connectivity.size())) {
    UInt conn = connectivity[i];
    for (UInt j : arange(nb_dofs)) {
      equations.push_back(conn * nb_dofs + j);
    }
  }
  
  for (UInt i : arange(total_nb_dofs)) {
    UInt row = equations[i];
    for (UInt j : arange(total_nb_dofs)) {
      UInt col = equations[j];
      global.add(row, col, local(i, j));
    }
  }
}
  
/* -------------------------------------------------------------------------- */
void Resolution::assembleStiffnessMatrix(GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  auto & stiffness =
      const_cast<SparseMatrix &>(model.getDOFManager().getMatrix("K"));

  auto & gaps = model.getGaps();
  auto & projections = model.getProjections();
  auto & normals = model.getNormals();
  
  UInt surface_dimension = spatial_dimension - 1;
  
  for (auto & element : model.getContactElements()) {

    auto nb_nodes  = element.getNbNodes();

    Real gap(gaps.begin()[element.slave]);
    Vector<Real> normal(normals.begin(spatial_dimension)[element.slave]);
    Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);
    
    Matrix<Real> covariant_basis(surface_dimension, spatial_dimension);
    GeometryUtils::covariantBasis(model.getMesh(), model.getContactDetector().getPositions(),
				  element.master, projection, covariant_basis);

    Vector<Real> delta_g(nb_nodes * spatial_dimension);
    ResolutionUtils::firstVariationNormalGap(element, projection, normal, delta_g);
       
    Matrix<Real> ddelta_g(nb_nodes * spatial_dimension, nb_nodes * spatial_dimension);
    ResolutionUtils::secondVariationNormalGap(element, covariant_basis,
					      projection, normal, gap, ddelta_g);
    
    Matrix<Real> k_n(nb_nodes * spatial_dimension, nb_nodes * spatial_dimension);
    computeNormalModuli(element, ddelta_g, delta_g, k_n);

    assembleLocalToGlobalMatrix(element, k_n, stiffness);
    
    /*Matrix<Real> m_alpha_beta(surface_dimension, surface_dimension);
    ResolutionUtils::computeMetricTensor(m_alpha_beta, element.tangents);

    // normal tangent moduli
        
    Array<Real> t_alpha(nb_nodes * spatial_dimension, surface_dimension);
    ResolutionUtils::computeTalpha(element, t_alpha);
    
    Array<Real> n_alpha(nb_nodes * spatial_dimension, surface_dimension);
    ResolutionUtils::computeNalpha(element, n_alpha);
    
    Array<Real> d_alpha(nb_nodes * spatial_dimension, surface_dimension);      
    ResolutionUtils::firstVariationNaturalCoordinate(element, d_alpha);

    computeNormalModuli(kc, n_alpha, d_alpha, delta_g, element);*/

    // frictional tangent moduli
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

    /*std::vector<UInt> equations;
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
      }*/
  }

  AKANTU_DEBUG_OUT();
}


} // namespace akantu
