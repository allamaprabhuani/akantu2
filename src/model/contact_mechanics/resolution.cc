/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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

    assembleLocalToGlobalArray(element, local_fn, model.getNormalForce());
    assembleLocalToGlobalArray(element, local_ft, model.getTangentialForce());
  }

  for (auto && [fc, ft, fn] :
       zip(make_view(model.getInternalForce(), spatial_dimension),
           make_view(model.getTangentialForce(), spatial_dimension),
           make_view(model.getNormalForce(), spatial_dimension))) {
    fc = ft + fn;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleLocalToGlobalArray(const ContactElement & element,
                                            Vector<Real> & local,
                                            Array<Real> & global) {

  auto get_connectivity = [&](auto & slave, auto & master) {
    const auto master_conn = model.getMesh().getConnectivity(master);
    Vector<Idx> elem_conn(master_conn.size() + 1);
    elem_conn[0] = slave;
    elem_conn.block(1, 0, master_conn.size(), 1) = master_conn;

    return elem_conn;
  };

  auto connectivity = get_connectivity(element.slave, element.master);

  auto nb_dofs = global.getNbComponent();
  auto nb_nodes = is_master_deformable ? connectivity.size() : 1;
  Real alpha = is_master_deformable ? 0.5 : 1.;

  MatrixProxy<Real> local_m(local.data(), nb_dofs, connectivity.size());
  auto global_it = make_view(global, nb_dofs).begin();

  for (Int i : arange(nb_nodes)) {
    global_it[connectivity(i)] = alpha * local_m(i);
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
    local_kn.zero();
    computeNormalModuli(element, local_kn);
    assembleLocalToGlobalMatrix(element, local_kn, global_stiffness);

    Matrix<Real> local_kt(nb_nodes * spatial_dimension,
                          nb_nodes * spatial_dimension);
    local_kt.zero();
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
    const auto master_conn = model.getMesh().getConnectivity(master);
    Vector<Idx> elem_conn(master_conn.size() + 1);
    elem_conn[0] = slave;
    elem_conn.block(1, 0, master_conn.size(), 1) = master_conn;

    return elem_conn;
  };

  auto connectivity = get_connectivity(element.slave, element.master);

  auto nb_dofs = spatial_dimension;
  UInt nb_nodes = is_master_deformable ? connectivity.size() : 1;
  UInt total_nb_dofs = nb_dofs * nb_nodes;

  std::vector<UInt> equations;
  for (Int i : arange(connectivity.size())) {
    UInt conn = connectivity[i];
    for (Int j : arange(nb_dofs)) {
      equations.push_back(conn * nb_dofs + j);
    }
  }

  for (Int i : arange(total_nb_dofs)) {
    UInt row = equations[i];
    for (Int j : arange(total_nb_dofs)) {
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
