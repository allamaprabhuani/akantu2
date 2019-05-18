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
void ResolutionPenalty::computeFrictionForce(Vector<Real> & force,
					     Array<Real>  & d_alpha,
					     Matrix<Real> & m_alpha_beta,
					     Vector<Real> & projection,
					     Real & gap) {

  auto traction = computeFrictionalTraction(m_alpha_beta, projection, gap);

  for (auto && values:
	 zip(traction,
	     make_view(d_alpha, d_alpha.size()))) {
    auto & t_s = std::get<0>(values);
    auto & d_s = std::get<1>(values);
    force += d_s * t_s;
  }
}

/* -------------------------------------------------------------------------- */
Vector<Real> ResolutionPenalty::computeFrictionalTraction(Matrix<Real>& m_alpha_beta,
							  Vector<Real>& xi,
							  Real & gap) {

  Real tn = gap * epsilon;
  tn = macaulay(tn);

  Vector<Real> previous_xi;
  
  auto delta_xi = xi - previous_xi;

  Vector<Real> previous_traction;
  Vector<Real> trial_traction;

  trial_traction.mul<false>(m_alpha_beta, delta_xi, epsilon);
  //trial_traction += previous_traction;
  
  auto trial_slip_function = trial_traction.norm() - mu * tn;

  Vector<Real> traction;
  if (trial_slip_function <= 0) {
    traction = trial_traction;
  }
  else{
    traction = mu * tn * trial_traction / trial_traction.norm(); 
  }  

  return traction;
}
  
  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentModuli(Matrix<Real> & kc, Vector<Real> & n,
					     Array<Real> & n_alpha, Array<Real> & d_alpha,
					     Matrix<Real> & m_alpha_beta,
					     Real & gap)  {
  computeNormalStiffness(kc, n, n_alpha, d_alpha, m_alpha_beta, gap);
  computeFrictionalStiffness(n, n_alpha, d_alpha, gap);
}

  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalStiffness(Matrix<Real> & ke, Vector<Real> & n,
					       Array<Real> & n_alpha, Array<Real> & d_alpha,
					       Matrix<Real> & /*m_alpha_beta*/, Real & gap) {

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
void ResolutionPenalty::computeFrictionalStiffness(Vector<Real> & /*n*/,
						   Array<Real>  & /*n_alpha*/,
						   Array<Real> & /*d_alpha*/,
						   Real & gap) {

  computeCommonModuli(gap);
  computeStickModuli();
  computeSlipModuli();

}


/* -------------------------------------------------------------------------- */
Matrix<Real> ResolutionPenalty::computeCommonModuli(Array<Real> & /*t_alpha_beta*/,
						    Array<Real> & /*n_alpha_beta*/,
						    Matrix<Real> & /*tangents*/
						    Vector<Real> & /* n */
						    Real & /*gap*/)  {

  Array<Real> kt_alpha(spatial_dimension -1, d_alpha.size() * d_alpha.size(),
			"k_T_alpha");
  for(auto && values :
	zip(tangents.transpose(),
	    make_view(kt_alpha, kt_alpha.size()),
	    make_view(t_alpha_beta, t_alpha_beta.size()),
	    make_view(n_alpha_beta, n_alpha_beta.size()))) {
    auto & tangent_s = std::get<0>(values);
    auto & kt_s      = std::get<1>(values);
    auto & t_alpha_s = std::get<2>(values);
    auto & n_alpha_s = std::get<3>(values);
    
    Matrix<Real> kt_s_mat(kt_s.storage(), d_alpha.size(), d_alpha.size())
    
    for(auto && tuple :
	  make_view(d_alpha, d_alpha.size())) {
      auto & d_s = std::get<0>(tuple);
      
      for(auto && entry :
	    make_view(d_alpha, d_alpha.size())) {
	auto & d_t = std::get<0>(entry);
	// compute constant
	Matrix<Real> tmp(d_t.size(), d_t.size());
	tmp.mul<false, true>(d_s, d_t);
	kt_s_mat += tmp;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ResolutionPenalty::computeStickModuli(Array<Real> & g_alpha,
						   Array<Real> & d_alpha,
						   Matrix<Real> & m_alpha_beta) {

  Matrix<Real> k_stick(d_alpha.size(), d_alpha.size());
  
  for (auto && values :
	 zip(arange(d_alpha.getNbComponent())
	     make_view(d_alpha, d_alpha.size()),
	     make_view(g_alpha, g_alpha.size()))) {
    auto & s   = std::get<0>(values); 
    auto & d_s = std::get<1>(values);
    auto & g_s = std::get<2>(values);

    Matrix<Real> ds_mat(d_s.storage(), d_s.size(), 1);
    Matrix<Real> gs_mat(g_s.storage(), g_s.size(), 1);

    Matrix<Real> tmp1(d_s.size(), d_s.size());
    tmp1.mul<false, true>(ds_mat, gs_mat);

    k_stick += tmp1;
    
    for(auto && tuple :
	  enumerate(make_view(d_alpha, d_alpha.size()))) {
      auto & t   = std::get<0>(tuple);
      auto & d_t = std::get<1>(tuple);
       
      Matrix<Real> dt_mat(d_t.storage(), d_t.size(), 1);
      Matrix<Real> tmp2(d_t.size(), d_t.size());

      tmp2.mul<false, true>(ds_mat, dt_mat);
      k_stick += tmp2 * m_alpha_beta(s, t);
    }

    k_stick *= epsilon_t;

    return k_stick;
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ResolutionPenalty::computeSlipModuli() {

}
  
INSTANTIATE_RESOLUTION(penalty, ResolutionPenalty);  
  
} // akantu
