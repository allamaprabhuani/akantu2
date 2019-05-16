/**
 * @file   resolution_penalty.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Jan 7 2019
 * @date last modification: Mon Jan 7 2019
 *
 * @brief  Specialization of the resolution class for the penalty method
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
#include "resolution_penalty.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */  
ResolutionPenalty::ResolutionPenalty(ContactMechanicsModel & model,
				     const ID & id)
    : Resolution(model, id) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

 
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::initialize() {
  this->registerParam("epsilon", epsilon, Real(0.), _pat_parsable | _pat_modifiable,
		      "Normal penalty parameter");
  this->registerParam("epsilon_t", epsilon_t, Real(0.), _pat_parsable | _pat_modifiable,
		      "Tangential penalty parameter");

}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalForce(Vector<Real> & force, Vector<Real> & n,
					   Real & gap) {

  force.clear();
  Real tn = gap * epsilon;
  tn = macaulay(tn);
  for (UInt i : arange(force.size())) {
      force[i] += tn * n[i];
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionForce(Vector<Real> & force, Array<Real> & d_alpha,
					     Real & /*gap*/) {

  Vector<Real> tractions(d_alpha.getNbComponent());
  /*computeFrictionalTraction(tractions);

  for (auto && values:
	 zip(tractions,
	     make_view(d_alpha, d_alpha.size()))) {
    auto & t_s = std::get<0>(values);
    auto & d_s = std::get<1>(values);
    force += d_s * t_s;
    }*/
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionTraction(Real & gap) {

  Real tn = gap * epsilon;
  tn = macaulay(tn);

  // delta_xi = xi- previous_xi
  //trial_traction = previous_traction + epsilon * metric_tensor * delta_xi;
  //trial_slip_function = trial_traction.mod - mu * tn;

  // Stick condition
  // if trial_slip_function <= 0 
  // traction = trial_traction;

  // Slip condition
  // else
  // traction = mu * tn * trial_traction / trial_traction.norm()'
  
  
}
  
  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentModuli(Matrix<Real> & kc, Vector<Real> & n,
					     Array<Real> & n_alpha, Array<Real> & d_alpha,
					     Matrix<Real> & surface_matrix,
					     Real & gap)  {
  computeNormalStiffness(kc, n, n_alpha, d_alpha, surface_matrix, gap);
  computeFrictionalStiffness(n, n_alpha, d_alpha, gap);
}

  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalStiffness(Matrix<Real> & ke, Vector<Real> & n,
					       Array<Real> & n_alpha, Array<Real> & d_alpha,
					       Matrix<Real> & /*surface_matrix*/, Real & gap) {

  Real tn = gap * epsilon;
  tn = macaulay(tn);

  Matrix<Real> n_mat(n.storage(), n.size(), 1);
  
  ke.mul<false, true>(n_mat, n_mat);
  ke *= epsilon * heaviside(gap);
  
  for (auto && values:
	 zip(make_view(n_alpha, n_alpha.size()),
	     make_view(d_alpha, d_alpha.size()))) {
    auto & n_s = std::get<0>(values);
    auto & d_s = std::get<1>(values);

    Matrix<Real> ns_mat(n_s.storage(), n_s.size(), 1);
    Matrix<Real> ds_mat(d_s.storage(), d_s.size(), 1);

    Matrix<Real> tmp1(n_s.size(), n_s.size());
    tmp1.mul<false, true>(ns_mat, ds_mat);

    Matrix<Real> tmp2(n_s.size(), n_s.size());
    tmp1.mul<false, true>(ds_mat, ns_mat);

    ke -= (tmp1 + tmp2) * tn;
  }

}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionalStiffness(Vector<Real> & n,
						   Array<Real>  & n_alpha, Array<Real> & d_alpha,
						   Real & gap) {

  computeCommonModuli();
  computeStickModuli();
  computeSlipModuli();

}


/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeCommonModuli(Real & /*gap*/)  {


}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeStickModuli() {

  // epsilon_t * 
  
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeSlipModuli() {

}
  
INSTANTIATE_RESOLUTION(penalty, ResolutionPenalty);  
  
} // akantu
