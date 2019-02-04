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
		      "Penalty parameter");
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalForce(Vector<Real> & force, Vector<Real> & n,
					   Real & gap) {
  Real tn = gap * epsilon;
  tn = macaulay(tn);
  for (UInt i : arange(force.size())) {
      force[i] += tn * n[i];
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionForce(Vector<Real> & force, Array<Real> & d_alpha,
					     Real & gap) {

  Vector<Real> tractions(d_alpha.getNbComponent());
  computeFrictionalTraction(tractions);

  for (auto && values:
	 zip(tractions,
	     make_view(d_alpha, d_alpha.size()))) {
    auto & t_s = std::get<0>(values);
    auto & d_s = std::get<1>(values);
    force += d_s * t_s;
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionalTraction(Vector<Real> & tractions) {

}
  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentModuli(Vector<Real> & n, Array<Real> & n_alpha,
					     Array<Real> & t_alpha, Array<Real> & d_alpha,
					     Matrix<Real> & surface_matrix, Real & gap)  {
  computeNormalStiffness(n, n_alpha, d_alpha, surface_matrix, gap);
  computeFrictionalStiffness(n, n_alpha, d_alpha, gap);
}

  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalStiffness(Vector<Real> & n,
					       Array<Real> & n_alpha, Array<Real> & d_alpha,
					       Matrix<Real> & surface_matrix, Real & gap) {

  Real tn = gap * epsilon;
  tn = macaulay(tn);

  Matrix<Real> ke(n_alpha.size(), n_alpha.size());
  Matrix<Real> n_mat(n.storage(), n.size(), 1);
  
  ke.mul<false, true>(n_mat, n_mat);
  ke *= epsilon * heaviside(gap);
  
  for (auto && entry1 :
	 enumerate(make_view(n_alpha, n_alpha.size()))) {
    
    auto & i = std::get<0>(entry1);
    auto & ni = std::get<1>(entry1);
    Matrix<Real> ni_mat(ni.storage(), ni.size(), 1);
    
    for (auto && entry2:
	   enumerate(make_view(n_alpha, n_alpha.size()))) {

      auto & j = std::get<0>(entry2);
      auto & nj = std::get<1>(entry2);
      Matrix<Real> nj_mat(nj.storage(), nj.size(), 1);

      Matrix<Real> tmp(nj.size(), nj.size());
      tmp.mul<false, true>(ni_mat, nj_mat);
      tmp *= surface_matrix[i, j];
     
      ke += tmp;
    }    
  }

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

    ke -= tmp1 + tmp2;
  }
 
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionalStiffness(Vector<Real> & n,
						   Array<Real>  & n_alpha, Array<Real> & d_alpha,
						   Real & gap) {
					       
}
   
  
INSTANTIATE_RESOLUTION(penalty, ResolutionPenalty);  
  
} // akantu
