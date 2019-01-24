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
void ResolutionPenalty::computeNormalForce(Vector<Real> & force, Array<Real> & n,
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
  ///needs to calculate tt

  Real tt = 0;
  for (auto && values:
	 zip(make_view(d_alpha, d_alpha.size()))) {
    auto & d =  std::get<0>(values);
    for (UInt i : arange(d.size())) {
      force[i] += d[i] * tt; 
    }
  }
      
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentModuli(Array<Real> & n, Array<Real> & n_alpha,
					     Array<Real> & t_alpha, Array<Real> & d_alpha,
					     Matrix<Real> & surface_matrix, Real & gap)  {
  computeNormalStiffness(n, n_alpha, d_alpha, gap);
  computeFrictionalStiffness(n, n_alpha, d_alpha, gap);
}

  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalStiffness(Array<Real> & n,
					       Array<Real>  & n_alpha, Array<Real> & d_alpha,
					       Real & gap) {
					       
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionalStiffness(Array<Real> & n,
						   Array<Real>  & n_alpha, Array<Real> & d_alpha,
						   Real & gap) {
					       
}
   
  
INSTANTIATE_RESOLUTION(penalty, ResolutionPenalty);  
  
} // akantu
