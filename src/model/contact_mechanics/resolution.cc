/**
 * @file   resolution.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jan 17 2019
 * @date last modification: Wed Apr 07 2021
 *
 * @brief  Implementation of common part of the contact resolution class
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
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
    : Parsable(ParserType::_contact_resolution, id), id(id),
      fem(model.getFEEngine()), model(model) {

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
                "Friction Coefficient");
  registerParam("is_master_deformable", is_master_deformable, bool(false),
                _pat_parsable | _pat_readable, "Is master surface deformable");
}

/* -------------------------------------------------------------------------- */
void Resolution::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);
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

  for (const auto & element : model.getContactElements()) {

    auto nb_nodes = element.getNbNodes();

    Vector<Real> local_fn(nb_nodes * spatial_dimension);
    computeNormalForce(element, local_fn);

    Vector<Real> local_ft(nb_nodes * spatial_dimension);
    computeTangentialForce(element, local_ft);

    Vector<Real> local_fc(nb_nodes * spatial_dimension);
    local_fc = local_fn + local_ft;

    assembleLocalToGlobalArray(element, local_fn, model.getNormalForce());
    assembleLocalToGlobalArray(element, local_ft, model.getTangentialForce());
    assembleLocalToGlobalArray(element, local_fc, model.getInternalForce());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleLocalToGlobalArray(const ContactElement & element,
                                            Vector<Real> & local,
                                            Array<Real> & global) {

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

  auto & surface_selector = model.getContactDetector().getSurfaceSelector();
  auto & slave_list = surface_selector.getSlaveList();
  auto & master_list = surface_selector.getMasterList();

  auto connectivity = get_connectivity(element.slave, element.master);

  UInt nb_dofs = global.getNbComponent();
  UInt nb_nodes = is_master_deformable ? connectivity.size() : 1;
  Real alpha = is_master_deformable ? 0.5 : 1.;

  for (UInt i : arange(nb_nodes)) {
    UInt n = connectivity[i];

    auto slave_result = std::find(slave_list.begin(), slave_list.end(), n);
    auto master_result = std::find(master_list.begin(), master_list.end(), n);

    for (UInt j : arange(nb_dofs)) {
      UInt offset_node = n * nb_dofs + j;
      global[offset_node] += alpha * local[i * nb_dofs + j];
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleStiffnessMatrix(GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  auto & global_stiffness =
      const_cast<SparseMatrix &>(model.getDOFManager().getMatrix("K"));

  for (const auto & element : model.getContactElements()) {

    auto nb_nodes = element.getNbNodes();

    Matrix<Real> local_kn(nb_nodes * spatial_dimension,
                          nb_nodes * spatial_dimension);
    computeNormalModuli(element, local_kn);
    assembleLocalToGlobalMatrix(element, local_kn, global_stiffness);

    Matrix<Real> local_kt(nb_nodes * spatial_dimension,
                          nb_nodes * spatial_dimension);
    computeTangentialModuli(element, local_kt);
    assembleLocalToGlobalMatrix(element, local_kt, global_stiffness);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleLocalToGlobalMatrix(const ContactElement & element,
                                             const Matrix<Real> & local,
                                             SparseMatrix & global) {

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

  auto nb_dofs = spatial_dimension;
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
void Resolution::beforeSolveStep() {}

/* -------------------------------------------------------------------------- */
void Resolution::afterSolveStep(__attribute__((unused)) bool converged) {}

} // namespace akantu
