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
  this->registerParam("epsilon_n", epsilon_n, Real(0.),
                      _pat_parsable | _pat_modifiable,
                      "Normal penalty parameter");
  this->registerParam("epsilon_t", epsilon_t, Real(0.),
                      _pat_parsable | _pat_modifiable,
                      "Tangential penalty parameter");
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalForce(const ContactElement & element,
					   Vector<Real> & force) {

  force.clear();

  auto & gaps = model.getGaps();
  auto & projections = model.getProjections();
  auto & normals = model.getNormals();

  auto surface_dimension = spatial_dimension - 1;
  
  Real gap(gaps.begin()[element.slave]);
  Vector<Real> normal(normals.begin(spatial_dimension)[element.slave]);
  Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);

  // compute normal traction
  Real p_n = macaulay(gap) * epsilon_n;

  // compute first variation of gap
  auto nb_nodes = element.getNbNodes();
  Vector<Real> delta_gap(nb_nodes * spatial_dimension);
  ResolutionUtils::firstVariationNormalGap(element, projection, normal, delta_gap);

  // compute normal force
  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());
  for (UInt i : arange(force.size())) 
    force[i] += delta_gap[i] * p_n  * nodal_area[element.slave];
  
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentialForce(const ContactElement & element,
					       Vector<Real> & force) {
  
  if (mu == 0)
    return;
  
  UInt surface_dimension = spatial_dimension - 1; 
  
  // compute tangents
  auto & projections = model.getProjections();
  Vector<Real> projection(projections.begin(surface_dimension)[element.slave]);
  
  auto & normals = model.getNormals();
  Vector<Real> normal(normals.begin(spatial_dimension)[element.slave]);
  
  auto & tangents = model.getTangents();
  Matrix<Real> covariant_basis(tangents.begin(surface_dimension,
					      spatial_dimension)[element.slave]);
  GeometryUtils::covariantBasis(model.getMesh(), model.getContactDetector().getPositions(),
				element.master, normal, projection, covariant_basis);
  
  // check for no-contact to contact condition
  auto & previous_master_elements = model.getPreviousMasterElements();
  auto & previous_element = previous_master_elements[element.slave];

  if (previous_element.type == _not_defined)
    return;
  
  // compute tangential traction using return map algorithm
  auto & tangential_tractions = model.getTangentialTractions();
  Vector<Real> tangential_traction(tangential_tractions.begin(surface_dimension)[element.slave]);
  this->computeTangentialTraction(element, covariant_basis,
				  tangential_traction);
  
  // compute first variation of natural coordinate  
  auto & gaps = model.getGaps();
  auto & gap = gaps.begin()[element.slave];

  auto nb_nodes  = element.getNbNodes();
  Array<Real> delta_xi(nb_nodes * spatial_dimension, surface_dimension);
  ResolutionUtils::firstVariationNaturalCoordinate(element, covariant_basis,
						   projection, normal, gap, delta_xi);

  // compute tangential force
  auto & nodal_area = const_cast<Array<Real> &>(model.getNodalArea());
  for (auto && values : zip(tangential_traction,
			    make_view(delta_xi, delta_xi.size()))) {
    auto & traction_alpha = std::get<0>(values);
    auto & delta_xi_alpha = std::get<1>(values);
    force += delta_xi_alpha * traction_alpha * nodal_area[element.slave];
  }
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentialTraction(const ContactElement & element,
					   const Matrix<Real> & covariant_basis,
					   Vector<Real> & traction_tangential) {

  UInt surface_dimension = spatial_dimension - 1;
  
  auto & gaps = model.getGaps();
  auto & gap = gaps.begin()[element.slave];

  // Return map algorithm is employed
  // compute trial traction
  Vector<Real> traction_trial(surface_dimension);
  this->computeTrialTangentialTraction(element, covariant_basis, traction_trial);

  // compute norm of trial traction
  Real traction_trial_norm = 0;
  auto inv_A = GeometryUtils::contravariantMetricTensor(covariant_basis);  
  for (auto i : arange(surface_dimension)) {
    for (auto j : arange(surface_dimension)) {
      traction_trial_norm += traction_trial[i] * traction_trial[j] * inv_A(i, j);
    }
  }
  traction_trial_norm = sqrt(traction_trial_norm);

  // check stick or slip condition
  auto & stick_or_slip = model.getStickSlip();
  auto & cond = stick_or_slip.begin()[element.slave];
  
  Real p_n = macaulay(gap) * epsilon_n;
  bool stick = (traction_trial_norm <= mu * p_n) ? true : false;

  if (stick) {
    cond = 1;
    computeStickTangentialTraction(element, traction_trial, traction_tangential);
  } else {
    cond = 2;
    computeSlipTangentialTraction(element, covariant_basis, traction_trial,
				  traction_tangential);
  }

}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTrialTangentialTraction(const ContactElement & element,
						       const Matrix<Real> & covariant_basis,
						       Vector<Real> & traction) {
  
  UInt surface_dimension = spatial_dimension - 1;
  
  /*auto & projections = model.getProjections();
  Vector<Real> contravariant_projection(projections.begin(surface_dimension)[element.slave]);

  auto & stick_projections = model.getStickProjections();
  Vector<Real> covariant_stick(stick_projections.begin(surface_dimension)[element.slave]);
  
  Matrix<Real> contravariant_basis(surface_dimension, spatial_dimension);
  GeometryUtils::contravariantBasis(covariant_basis, contravariant_basis);

  Vector<Real> covariant_projection(surface_dimension);
  for (auto && values : zip(covariant_projection,
			    contravariant_projection,
			    contravariant_basis.transpose())) {
    auto & temp = std::get<0>(values);
    Vector<Real> contravariant_alpha(std::get<2>(values));

    temp = contravariant_alpha.dot(contravariant_alpha);
    temp *= std::get<1>(values);
  }

  auto covariant_slip = covariant_projection - covariant_stick;
  traction.mul<true>(contravariant_basis, covariant_slip, epsilon_t);*/

  auto & projections =  model.getProjections();
  Vector<Real> current_projection(projections.begin(surface_dimension)[element.slave]);

  auto & previous_projections = model.getPreviousProjections();
  Vector<Real> previous_projection(previous_projections.begin(surface_dimension)[element.slave]);

  // method from Laursen et. al.
  /*auto covariant_metric_tensor = GeometryUtils::covariantMetricTensor(covariant_basis);
  auto increment_projection = current_projection - previous_projection;

  traction.mul<false>(covariant_metric_tensor, increment_projection, epsilon_t);
  
  auto & previous_tangential_tractions = model.getPreviousTangentialTractions();
  Vector<Real> previous_traction(previous_tangential_tractions.begin(surface_dimension)[element.slave]);
  traction = previous_traction + traction;*/

  // method from Schweizerhof
  auto covariant_metric_tensor = GeometryUtils::covariantMetricTensor(covariant_basis);

  auto & previous_tangential_tractions = model.getPreviousTangentialTractions();
  Vector<Real> previous_traction(previous_tangential_tractions.begin(surface_dimension)[element.slave]);

  auto & previous_tangents = model.getPreviousTangents();
  Matrix<Real> previous_covariant_basis(previous_tangents.begin(surface_dimension,
								spatial_dimension)[element.slave]);
  auto previous_contravariant_metric_tensor =
    GeometryUtils::contravariantMetricTensor(previous_covariant_basis);

  auto current_tangent = covariant_basis.transpose();
  auto previous_tangent = previous_covariant_basis.transpose();
  
  for (auto alpha :arange(surface_dimension)) {
    Vector<Real> tangent_alpha(current_tangent(alpha));
    for (auto gamma : arange(surface_dimension)) {
      for (auto beta : arange(surface_dimension)) {
	Vector<Real> tangent_beta(previous_tangent(beta));
	auto t_alpha_t_beta = tangent_beta.dot(tangent_alpha);
	traction[alpha] += previous_traction[gamma]*previous_contravariant_metric_tensor(gamma, beta)*t_alpha_t_beta;
      }
    }
  }
  
  auto & previous_master_elements = model.getPreviousMasterElements();
  auto & previous_element = previous_master_elements[element.slave];

  Vector<Real> previous_real_projection(spatial_dimension);
  GeometryUtils::realProjection(model.getMesh(), model.getContactDetector().getPositions(),
				previous_element, previous_projection,
				previous_real_projection);
      
  Vector<Real> current_real_projection(spatial_dimension);
  GeometryUtils::realProjection(model.getMesh(), model.getContactDetector().getPositions(),
				element.master, current_projection,
				current_real_projection);

  auto increment_real = current_real_projection - previous_real_projection;
  Vector<Real> increment_xi(surface_dimension);

  auto contravariant_metric_tensor = GeometryUtils::contravariantMetricTensor(covariant_basis);

  // increment in natural coordinate
  for (auto beta : arange(surface_dimension)) {
    for (auto gamma: arange(surface_dimension)) {
      auto temp = increment_real.dot(current_tangent(gamma));
      temp *= contravariant_metric_tensor(beta, gamma);
      increment_xi[beta] += temp;
    }
  }

  Vector<Real> temp(surface_dimension);
  temp.mul<false>(covariant_metric_tensor, increment_xi, epsilon_t);

  traction -= temp;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeStickTangentialTraction(const ContactElement & /*element*/,
						       Vector<Real> & traction_trial,
						       Vector<Real> & traction_tangential) {
  traction_tangential = traction_trial;
}
  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeSlipTangentialTraction(const ContactElement & element,
						      const Matrix<Real> & covariant_basis,
						      Vector<Real> & traction_trial,
						      Vector<Real> & traction_tangential) {
  UInt surface_dimension = spatial_dimension - 1; 

  auto & gaps = model.getGaps();
  auto & gap = gaps.begin()[element.slave];

  // compute norm of trial traction
  Real traction_trial_norm = 0;
  auto inv_A = GeometryUtils::contravariantMetricTensor(covariant_basis); 
  for (auto i : arange(surface_dimension)) {
    for (auto j : arange(surface_dimension)) {
      traction_trial_norm += traction_trial[i] * traction_trial[j] * inv_A(i, j);
    }
  }
  traction_trial_norm = sqrt(traction_trial_norm);
  
  auto slip_direction = traction_trial;
  slip_direction /= traction_trial_norm;

  Real p_n = epsilon_n * macaulay(gap);
  traction_tangential = slip_direction;
  traction_tangential *= p_n;

  /*
  auto & stick_projections = model.getStickProjections();
  Vector<Real> covariant_stick(stick_projections.begin(surface_dimension)[element.slave]);

  
  auto slip = macaulay(traction_trial_norm - mu * p_n);
  slip /= epsilon_t;

  Vector<Real> covariant_slip(surface_dimension);
  for (auto && values : zip(covariant_slip,
			    covariant_basis.transpose())) {

    auto & slip_alpha = std::get<0>(values);
    Vector<Real> covariant_alpha(std::get<1>(values));
    slip_alpha = slip_direction.dot(covariant_alpha);
  }
  
  covariant_stick += slip * covariant_slip;*/
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeNormalModuli(const ContactElement & element,
					    Matrix<Real> & ddelta_g,
					    Vector<Real> & delta_g, Matrix<Real> & stiffness) {

  auto & gaps = model.getGaps();
  auto & gap = gaps.begin()[element.slave];

  auto & nodal_areas = model.getNodalArea();
  auto & nodal_area = nodal_areas.begin()[element.slave];
  
  Matrix<Real> tmp(delta_g.storage(), delta_g.size(), 1);
  Matrix<Real> mat_delta_g(delta_g.size(), delta_g.size());

  Real heaviside_gap = heaviside(gap);
  mat_delta_g.mul<false, true>(tmp, tmp, heaviside_gap);

  Real macaulay_gap = macaulay(gap);
  ddelta_g *= macaulay_gap;
  
  stiffness += mat_delta_g + ddelta_g;
  stiffness *= epsilon_n * nodal_area;
}
  
/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeTangentialModuli(const ContactElement & /*element*/,
						Matrix<Real> & /*ddelta_g*/,
						Vector<Real> & /*delta_g*/,
						Matrix<Real> & /*stiffness*/){

  if (mu == 0) 
    return;
  
}
  
/* -------------------------------------------------------------------------- */
/*void ResolutionPenalty::computeNormalModuli(Matrix<Real> & ke,
                                            Array<Real> & n_alpha,
                                            Array<Real> & d_alpha,
                                            Vector<Real> & n,
                                            ContactElement & element) {

  auto & gaps = model.getGaps();
  auto & gap = gaps.begin()[element.slave];
  
  Real tn = gap * epsilon_n;
  tn = macaulay(tn);

  Matrix<Real> n_mat(n.storage(), n.size(), 1);
  ke.mul<false, true>(n_mat, n_mat);
  ke *= epsilon_n * heaviside(gap);
  
  for (auto && values : zip(make_view(n_alpha, n_alpha.size()),
                            make_view(d_alpha, d_alpha.size()))) {
    auto & n_s = std::get<0>(values);
    auto & d_s = std::get<1>(values);

    Matrix<Real> ns_mat(n_s.storage(), n_s.size(), 1);
    Matrix<Real> ds_mat(d_s.storage(), d_s.size(), 1);

    Matrix<Real> tmp1(n_s.size(), n_s.size());
    tmp1.mul<false, true>(ns_mat, ds_mat);

    Matrix<Real> tmp2(n_s.size(), n_s.size());
    tmp2.mul<false, true>(ds_mat, ns_mat);

    ke -= (tmp1 + tmp2) * tn;
  }

  auto surface_dimension = spatial_dimension - 1;
  Matrix<Real> m_alpha_beta(surface_dimension, surface_dimension);
  ResolutionUtils::computeMetricTensor(m_alpha_beta, element.tangents);

  for (auto && values :
       zip(arange(surface_dimension), make_view(n_alpha, n_alpha.size()))) {

    auto & s = std::get<0>(values);
    auto & n_s = std::get<1>(values);

    Matrix<Real> ns_mat(n_s.storage(), n_s.size(), 1);

    for (auto && tuple :
         zip(arange(surface_dimension), make_view(n_alpha, n_alpha.size()))) {
      auto & t = std::get<0>(tuple);
      auto & n_t = std::get<1>(tuple);

      Matrix<Real> nt_mat(n_t.storage(), n_t.size(), 1);

      Matrix<Real> tmp3(n_s.size(), n_s.size());
      tmp3.mul<false, true>(ns_mat, nt_mat);

      tmp3 *= m_alpha_beta(s, t) * tn * element.gap;
      ke += tmp3;
    }
  }
}*/

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::computeFrictionalModuli(
    Matrix<Real> & /*ke*/, Array<Real> & t_alpha_beta,
    Array<Real> & n_alpha_beta, Array<Real> & /*n_alpha*/,
    Array<Real> & d_alpha, Matrix<Real> & phi, Vector<Real> & n,
    ContactElement & element) {

  /*auto k_common =
      computeCommonModuli(t_alpha_beta, n_alpha_beta, d_alpha, n, element);

  const auto & type = element.master.type;
  const auto & conn = element.connectivity;

  auto surface_dimension = Mesh::getSpatialDimension(type);
  auto spatial_dimension = surface_dimension + 1;

  Matrix<Real> m_alpha_beta(surface_dimension, surface_dimension);
  ResolutionUtils::computeMetricTensor(m_alpha_beta, element.tangents);

  Array<Real> g_alpha(conn.size() * spatial_dimension, surface_dimension);
  ResolutionUtils::computeGalpha(g_alpha, t_alpha_beta, d_alpha, phi, element);*/

  /*Matrix<Real> k_t;
  bool stick = computeFrictionalTraction(m_alpha_beta, element);

  if(stick)
    k_t = computeStickModuli(g_alpha, d_alpha, m_alpha_beta);
  else
    k_t = computeSlipModuli(g_alpha, d_alpha, m_alpha_beta, element);*/
}

/* -------------------------------------------------------------------------- */
Array<Real> ResolutionPenalty::computeCommonModuli(
    Array<Real> & /*t_alpha_beta*/, Array<Real> & /*n_alpha_beta*/,
    Array<Real> & d_alpha, Vector<Real> & /*n*/, ContactElement & /*element*/) {

  Array<Real> kt_alpha(spatial_dimension - 1, d_alpha.size() * d_alpha.size(),
                       "k_T_alpha");

  // auto t_alpha_beta_size = t_alpha_beta.size() * (spatial_dimension - 1);

  // auto & tangents = element.tangents;
  /*for(auto && values :
    zip(tangents.transpose(),
        make_view(kt_alpha, kt_alpha.size()),
        make_view(t_alpha_beta, t_alpha_beta_size),
        make_view(n_alpha_beta, n_alpha_beta_size))) {
    auto & tangent_s = std::get<0>(values);
    auto & kt_s      = std::get<1>(values);
    auto & t_alpha_s = std::get<2>(values);
    auto & n_alpha_s = std::get<3>(values);

    Matrix<Real> kt_s_mat(kt_s.storage(), d_alpha.size(), d_alpha.size());

    // loop over beta
    for(auto && tuple :
      zip(make_view(d_alpha, d_alpha.size()),
          make_view(n_alpha_ ))) {
      auto & d_s = std::get<0>(tuple);
      Matrix<Real> tmp(d_s.size(), d_s.size());

      // loop over gamma
      for(auto && entry :
        make_view(d_alpha, d_alpha.size())) {
    auto & d_t = std::get<0>(entry);


    // compute constant
    Matrix<Real> tmp2(d_t.size(), d_t.size());
    tmp2.mul<false, true>(d_s, d_t);
    kt_s_mat += tmp2;
      }
      }
      }*/

  return kt_alpha;
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ResolutionPenalty::computeStickModuli(
    Array<Real> & g_alpha, Array<Real> & d_alpha, Matrix<Real> & m_alpha_beta) {

  /*Matrix<Real> k_stick(d_alpha.size(), d_alpha.size());

  for (auto && values :
       zip(arange(d_alpha.getNbComponent()), make_view(d_alpha, d_alpha.size()),
           make_view(g_alpha, g_alpha.size()))) {
    auto & s = std::get<0>(values);
    auto & d_s = std::get<1>(values);
    auto & g_s = std::get<2>(values);

    Matrix<Real> ds_mat(d_s.storage(), d_s.size(), 1);
    Matrix<Real> gs_mat(g_s.storage(), g_s.size(), 1);

    Matrix<Real> tmp1(d_s.size(), d_s.size());
    tmp1.mul<false, true>(ds_mat, gs_mat);

    k_stick += tmp1;

    for (auto && tuple : enumerate(make_view(d_alpha, d_alpha.size()))) {
      auto & t = std::get<0>(tuple);
      auto & d_t = std::get<1>(tuple);

      Matrix<Real> dt_mat(d_t.storage(), d_t.size(), 1);
      Matrix<Real> tmp2(d_t.size(), d_t.size());

      tmp2.mul<false, true>(ds_mat, dt_mat);
      k_stick += tmp2 * m_alpha_beta(s, t);
    }
  }

  k_stick *= epsilon_t;
  return k_stick;*/
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ResolutionPenalty::computeSlipModuli(
						  Array<Real> & /*g_alpha*/, Array<Real> & d_alpha,
						  Matrix<Real> & /*m_alpha_beta*/, ContactElement & /*element*/) {

  /*Real tn = element.gap * epsilon;
  tn = macaulay(tn);

  Real factor;
  factor = epsilon_t * mu * tn;

  auto p_t = element.traction;
  p_t /= p_t.norm();*/

  Matrix<Real> k_slip(d_alpha.size(), d_alpha.size());

  /*
  // loop for alpha
  for(auto && value :
    make_view(d_alpha, d_alpha.size())) {
    auto & d_s = std::get<0>(value);

    // loop for beta
    for(auto && tuple :
      zip(arange(spatial_dimension - 1),
          make_view(d_alpha, d_alpha.size()),
          make_view(g_alpha, g_alpha.size()))) {
      auto & beta = std::get<0>(tuple);
      auto & d_beta = std::get<1>(tuple);
      auto & g_beta = std::get<2>(tuple);

      // loop for gamma
      for(auto && entry :
        zip(arange(spatial_dimension - 1),
        make_view(d_alpha, d_alpha.size()))) {
    auto & gamma = std::get<0>(entry);
    auto & d_gamma = std::get<1>(entry);

      }
    }


    }*/

  return k_slip;
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::beforeSolveStep() {
}

/* -------------------------------------------------------------------------- */
void ResolutionPenalty::afterSolveStep(bool converged) {

  auto method = model.getAnalysisMethod();
  if (method == _explicit_lumped_mass) {
    return ;
  }
  
  auto & K =
      const_cast<SparseMatrix &>(model.getDOFManager().getMatrix("K"));

  auto k_min = K.min();
  auto roundoff_error = 1e-17;

  const auto blocked_dofs = model.getDOFManager().getBlockedDOFs("displacement");
  Real nb_unknowns = 0;
  for (auto & bld : make_view(blocked_dofs)) {
    if (not bld)
      nb_unknowns++;
  }

  auto max_epsilon_n = k_min / sqrt(nb_unknowns * roundoff_error);
  if (epsilon_n > max_epsilon_n)
    epsilon_n = max_epsilon_n;
  
}
  
INSTANTIATE_RESOLUTION(penalty, ResolutionPenalty);

} // namespace akantu
