/**
 * @file   newmark-beta.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 23 2015
 * @date last modification: Wed Mar 27 2019
 *
 * @brief  implementation of the  newmark-@f$\beta@f$ integration  scheme.  This
 * implementation is taken from Méthodes  numériques en mécanique des solides by
 * Alain Curnier \note{ISBN: 2-88074-247-1}
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "newmark-beta-NL-beam.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"
// #include "angle_tool.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NewmarkBetaNLBeam::NewmarkBetaNLBeam(DOFManager & dof_manager, const ID & dof_id)
    : CentralDifference(dof_manager, dof_id) {
}

/* -------------------------------------------------------------------------- */

void NewmarkBetaNLBeam::predictor(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                            Array<Real> & u_dot_dot,
                            const Array<bool> & blocked_dofs) const {
  
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.size();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real * u_dot_val = u_dot.storage();
  Real * u_dot_dot_val = u_dot_dot.storage();
  bool * blocked_dofs_val = blocked_dofs.storage();


  for (UInt d = 0; d < nb_degree_of_freedom; d++) {
    if (!(*blocked_dofs_val)) {

      *u_dot_val = *u_dot_val + delta_t / 2 * *u_dot_dot_val;
    }
    u_dot_val++;
    u_dot_dot_val++;
    blocked_dofs_val++;
  }
  UInt dim = u.getNbComponent() / 2;
  for(auto && data : zip(make_view(blocked_dofs, dim, 2), make_view(u,dim, 2), make_view(u_dot,dim,2))) {
    Vector<bool> blocked_u = std::get<0>(data)(0);
    Vector<bool> blocked_theta = std::get<0>(data)(1);
    Vector<Real> disp = std::get<1>(data)(0);
    Vector<Real> theta = std::get<1>(data)(1);
    Vector<Real> vel = std::get<2>(data)(0);
    Vector<Real> rv = std::get<2>(data)(1);

    disp += delta_t * vel;
    //theta = log_map(exp_map(delta_t * this->rv[nd]) * exp_map(this->theta[nd]));
}

  AKANTU_DEBUG_OUT();
  
}


/* -------------------------------------------------------------------------- */
void NewmarkBetaNLBeam::corrector(Real delta_t,
                            Array<Real> & u, Array<Real> & u_dot,
                            Array<Real> & u_dot_dot,
                            const Array<bool> & blocked_dofs) const {
  
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.size();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real * u_dot_val = u_dot.storage();
  Real * u_dot_dot_val = u_dot_dot.storage();
  bool * blocked_dofs_val = blocked_dofs.storage();


  for (UInt d = 0; d < nb_degree_of_freedom; d++) {
    if (!(*blocked_dofs_val)) {

      *u_dot_val = *u_dot_val + delta_t / 2 * *u_dot_dot_val;
    }
    u_dot_val++;
    u_dot_dot_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
  
}

} // namespace akantu
