/**
 * @file   contact_manager_implicit.cc
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @brief  Implicit contact implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "implicit_contact_manager.hh"

__BEGIN_AKANTU__


template <>
bool ContactData <2, SolidMechanicsModel>::computeTangentAndResidual(real_map &lambda_new, search_type *fn, Real& error) {
	constexpr UInt dim = 2;
  
	const Array <Real> & position = model_.getCurrentPosition();
	const Real tol = (*this)[Uzawa_tol];
  
	// get global stiffness matrix and force vector
	SparseMatrix        &K     = model_.getStiffnessMatrix();
	Array <Real>        &F     = model_.getResidual();
	const Array <Int>   &eqnum = model_.getDOFSynchronizer().getLocalDOFEquationNumbers();
  
	static bool auto_flag = true;
	if (auto_flag) {
		auto_flag = false;
		if (!(*this)[Automatic_penalty_parameter]) {
			for (auto it = sm_.begin(); it != sm_.end(); ++it)
        penalty_[it->first] = (*this)[Epsilon];
		}
		// else get penalty values automatically
		else
      getPenaltyValues();
	}
  
	Real lm_diff = 0;
	Real lm_max = 0;
  
	auto it = sm_.begin();
	while (it != sm_.end()) {
		auto slave = it->first;
    
		Real epsilon = (*this)[Alpha] * penalty_[slave];
    
		// get slave point
		point_type s(&position(slave));
    
		auto master = it->second;
		bool no_master = master == element_type();
    
		// if node lies outside triangle
		if (no_master || !has_projection(s,point_type(&position(master.node(0))),
                                             point_type(&position(master.node(1))))) {
      
			auto r = fn->search(&position(slave));
      
			// try to find a new master
			if (r != -1) {
				it->second = master = element_type(model_, _segment_2, r);
			}
			// else remove master-slave pair from simulation
			else {

				master = element_type();
        
				gaps_.erase(slave);
				lambda_new.erase(slave);
				++it;
				continue;
			}
		}
    
		assert(master.type == _segment_2);
    
		// compute the point on the surface
		point_type a(&position(master.node(0)));
		point_type b(&position(master.node(1)));

    //      point_type p = closest_point_to_segment(s,a,b);
    //      vector_type xi = invert_map<vector_type, point_type>(p, a, b);

		Distance_minimizator <2, _segment_2> dm(s, master.coordinates());
		//    Distance_minimzer<_segment_2> dm(s, master.coordinates());
		vector_type xi(1, dm.master_coordinates()[0]);
		point_type p = dm.point();
    
		// compute normal
		vector_type nu = master.normal();
		point_type nup(static_cast <const Real *>(nu.data()));
    
		// compute and save gap
		Real gap = -(nup * (s - p));
		gaps_[slave] = gap;
    
		Real lambda_hat = multipliers_[slave] + epsilon * gap;
    
		if (lambda_hat < 0) {
			// increase iterator
			++it;
			// save value of lambda
			lambda_new[slave] = 0;
			continue;
		}
    
		Real s1 = epsilon * Heaviside(lambda_hat);
		Real s2 = Macauley(lambda_hat); // max(0,lambda_hat)
    
		// NEVER EXIT IF GAP IS NEGATIVE
		//    if (gap < 0)
		//      continue;
    
		std::vector <UInt> conn(master.numNodes() + 1); // 1 slave (not hardcoded)
		conn[0] = slave;
		for (UInt i = 0; i < master.numNodes(); ++i)
      conn[1 + i] = master.node(i);
    
		// evaluate shape functions at slave master coordinate
		vector_type sh(master.numNodes());
		compute_shapes(xi, sh);

		// compute vector N
		vector_type N(dim * (master.numNodes() + 1));
		for (UInt i = 0; i < dim; ++i) {
			N[i] = nu[i];
			for (UInt j = 0; j < master.numNodes(); ++j)
        N[(1 + j) * dim + i] = -nu[i] * sh[j];
		}

		// compute vector T
		vector_type DN(sh.size());
		compute_shape_derivatives(xi, DN);
		point_type tau = DN[0] * a + DN[1] * b;

		vector_type T(dim * (master.numNodes() + 1));
		for (UInt i = 0; i < dim; ++i) {
			T[i] = tau[i];
			for (UInt j = 0; j < master.numNodes(); ++j)
      T[(1 + j) * dim + i] = -tau[i] * sh[j];
		}

		// compute N1
		vector_type N1(dim * (master.numNodes() + 1));
		for (UInt i = 0; i < dim; ++i) {
			for (UInt j = 0; j < master.numNodes(); ++j)
      N1[(1 + j) * dim + i] = -nu[i] * DN[j];
		}

		// compute m11
		Real m11 = tau * tau;

		// compute D1
		vector_type D1 = T + gap * N1;
		D1 *= 1. / m11;

		// Note: N1bar = N1 - k11*D1, but since k11 = 0 for 2D, then
		// N1bar = N1
		vector_type &N1bar = N1;

		// stiffness matrix (only non-zero terms for 2D implementation)
		matrix_type kc = s1 * N * transpose(N);      // first term
		kc += (s2 * gap * m11) * N1bar * transpose(N1bar); // second term
		kc -= s2 * D1 * transpose(N1);               // sixth term
		kc -= s2 * N1 * transpose(D1);               // eight term

		// residual vector
		vector_type fc = s2 * N;

		assert(kc.rows() == fc.size());

		// assemble local components into global matrix and vector
		std::vector <UInt> eq;
		for (UInt i = 0; i < conn.size(); ++i)
      for (UInt j = 0; j < dim; ++j)
        eq.push_back(eqnum(conn[i] * dim + j));
    
		for (UInt i = 0; i < kc.rows(); ++i) {
			F[eq[i]] += fc(i);
			for (UInt j = i; j < kc.columns(); ++j) {
				K.addToProfile(eq[i], eq[j]);
				K(eq[i], eq[j]) += kc(i, j);
			}
		}
    
		// update multiplier
		lambda_new[slave] = s2;
    
		Real lm_old = multipliers_[slave];
		lm_max += lm_old * lm_old;
		lm_old -= s2;
		lm_diff += lm_old * lm_old;
    
		// increase iterator
		++it;
	}
  
	if (lm_max < tol) {
		error = sqrt(lm_diff);
		return sqrt(lm_diff) < tol;
	}
  
	error = sqrt(lm_diff / lm_max);
	return sqrt(lm_diff / lm_max) < tol;
}



template <>
bool ContactData <3, SolidMechanicsModel>::computeTangentAndResidual(real_map &lambda_new, search_type *fn, Real& error) {
	constexpr UInt dim = 3;
  
	const Array <Real> & position = model_.getCurrentPosition();
	const Real tol = (*this)[Uzawa_tol];
  
	// get global stiffness matrix and force vector
	SparseMatrix        &K     = model_.getStiffnessMatrix();
	Array <Real>        &F     = model_.getResidual();
	const Array <Int>   &eqnum = model_.getDOFSynchronizer().getLocalDOFEquationNumbers();
  
	static bool auto_flag = true;
	if (auto_flag) {
		auto_flag = false;
		if (!(*this)[Automatic_penalty_parameter]) {
			for (auto it = sm_.begin(); it != sm_.end(); ++it)
      penalty_[it->first] = (*this)[Epsilon];
		}
		// else get penalty values automatically
		else
    getPenaltyValues();
	}
  
  
	//    K.saveMatrix("something.mtx");
  
	Real lm_diff = 0;
	Real lm_max = 0;
  
	auto it = sm_.begin();
	while (it != sm_.end()) {
		auto slave = it->first;
    
		Real epsilon = (*this)[Alpha] * penalty_[slave];
    
		// get slave point
		point_type s(&position(slave));
    
		auto master = it->second;
		bool no_master = master == element_type();
    
		// if node lies outside triangle
		if (no_master || !point_has_projection_to_triangle(s,
		                                                   point_type(&position(master.node(0))),
		                                                   point_type(&position(master.node(1))),
		                                                   point_type(&position(master.node(2))))) {
      
			auto r = fn->search(&position(slave));

			// try to find a new master
			if (r != -1) {
				it->second = master = element_type(model_, _triangle_3, r);
			}
			// else remove master-slave pair from simulation
			else {
        //				if (this->niter_ > 10)
        //					cout << "- Info: default master for slave "  << endl;
        
				master = element_type();
        
				gaps_.erase(slave);
				lambda_new.erase(slave);
        
        //        auto eit = it;
				++it;
        //        sm_.erase(eit);
				continue;
			}
		}
    
		assert(master.type == _triangle_3);
    
		// compute the point on the surface
		point_type a(&position(master.node(0)));
		point_type b(&position(master.node(1)));
		point_type c(&position(master.node(2)));
    
    
		//    point_type p = closest_point_to_triangle(s, a, b, c);
		//    vector_type xi = invert_map<vector_type, point_type>(p, a, b, c);
    
		Distance_minimizator <3, _triangle_3> dm(s, master.coordinates());
		vector_type xi(2);
		xi[0] = dm.master_coordinates()[0];
		xi[1] = dm.master_coordinates()[1];
		point_type p = dm.point();
    
		// compute normal
		vector_type nu = master.normal();
		point_type nup(static_cast <const Real *>(nu.data()));
    
		// compute and save gap
		Real gap = -(nup * (s - p));
		gaps_[slave] = gap;
    
		Real lambda_hat = multipliers_[slave] + epsilon * gap;
    
		if (lambda_hat < 0) {
			// increase iterator
			++it;
			// save value of lambda
			lambda_new[slave] = 0;
			continue;
		}
    
		Real s1 = epsilon * Heaviside(lambda_hat);
		Real s2 = Macauley(lambda_hat); // max(0,lambda_hat)
		Real s3 = s2 * gap;
    
		// NEVER EXIT IF GAP IS NEGATIVE
		//    if (gap < 0)
		//      continue;
    
		std::vector <UInt> conn(master.numNodes() + 1); // 1 slave (not hardcoded)
		conn[0] = slave;
		for (UInt i = 0; i < master.numNodes(); ++i)
    conn[1 + i] = master.node(i);
    
    
		// evaluate shape functions at slave master coordinate
		vector_type sh(master.numNodes());
		InterpolationElement <_itp_lagrange_triangle_3>::computeShapes(xi, sh);
    
		// compute vector N
		vector_type N(dim * (master.numNodes() + 1));
		for (UInt i = 0; i < dim; ++i) {
			N[i] = nu[i];
			for (UInt j = 0; j < master.numNodes(); ++j)
      N[(1 + j) * dim + i] = -nu[i] * sh[j];
		}
    
		matrix_type DN(2, master.numNodes());
		InterpolationElement <_itp_lagrange_triangle_3>::computeDNDS(xi, DN);
    
		point_type tau1 = DN(0, 0) * a + DN(0, 1) * b + DN(0, 2) * c;
		point_type tau2 = DN(1, 0) * a + DN(1, 1) * b + DN(1, 2) * c;
    
		vector_type nucheck(3);
    
		Math::vectorProduct3(&tau1[0], &tau2[0], &nucheck[0]);
		Math::normalize3(&nucheck[0]);
    
		if ((nucheck - nu)().norm() > 1.0e-10) {
			cout << "*** ERROR *** Normal failed" << endl;
			cout << "nu1: " << nu << endl;
			cout << "nu2: " << nucheck << endl;
			exit(1);
		}
    
		// compute vectors T1, T2, N1, N2
		size_t vsize = dim * (master.numNodes() + 1);
		vector_type T1(vsize), T2(vsize), N1(vsize), N2(vsize);
		for (UInt i = 0; i < dim; ++i) {
			T1[i] = tau1[i];
			T2[i] = tau2[i];
			for (UInt j = 0; j < master.numNodes(); ++j) {
				T1[(1 + j) * dim + i] = -tau1[i] * sh[j];
				T2[(1 + j) * dim + i] = -tau2[i] * sh[j];
				N1[(1 + j) * dim + i] = -nu[i] * DN(0u, j);
				N2[(1 + j) * dim + i] = -nu[i] * DN(1u, j);
			}
		}
    
		// compute matrix A = m + k*g  (but kappa is zero for linear elements)
		matrix_type A(2, 2);
		A(0, 0) = tau1 * tau1;
		A(1, 1) = tau2 * tau2;
		A(0, 1) = tau1 * tau2;
    
		Real detA = A(0, 0) * A(1, 1) - A(0, 1) * A(0, 1);
    
		// compute vectors D1, D2
		vector_type D1 = (1 / detA) * (A(1, 1) * (T1 + gap * N1)() - A(0, 1) * (T2 + gap * N2)())();
		vector_type D2 = (1 / detA) * (A(0, 0) * (T2 + gap * N2)() - A(0, 1) * (T1 + gap * N1)())();
    
		// Note: N1bar = N1 - k11*D1, but since k11 = 0 for linear elements, then
		// N1bar = N1, N2bar = N2
		vector_type &N1bar = N1;
		vector_type &N2bar = N2;
    
		// stiffness matrix (only non-zero terms for 3D implementation with linear elements)
		matrix_type kc = s1 * N * transpose(N);   // first term
		kc += (s3 * A(0, 0)) * N1bar * transpose(N1bar); // second term
		matrix_type tmp = N1bar * transpose(N2bar);
		tmp += N2bar * transpose(N1bar);
		kc += (s3 * A(0, 1)) * tmp;               // third and fourth terms
		kc += (s3 * A(1, 1)) * N2bar * transpose(N2bar); // fifth term
		kc -= s2 * D1 * transpose(N1);            // sixth term
		kc -= s2 * D2 * transpose(N2);            // seventh term
		kc -= s2 * N1 * transpose(D1);            // eight term
		kc -= s2 * N2 * transpose(D2);            // ninth term
    
		// residual vector
		vector_type fc = s2 * N;
    
		assert(kc.rows() == fc.size());
    
		// assemble local components into global matrix and vector
		std::vector <UInt> eq;
		for (UInt i = 0; i < conn.size(); ++i)
    for (UInt j = 0; j < dim; ++j)
    eq.push_back(eqnum(conn[i] * dim + j));
    
		for (UInt i = 0; i < kc.rows(); ++i) {
			F[eq[i]] += fc(i);
			for (UInt j = i; j < kc.columns(); ++j) {
				K.addToProfile(eq[i], eq[j]);
				K(eq[i], eq[j]) += kc(i, j);
			}
		}
    
		// update multiplier
		lambda_new[slave] = s2;
    
		Real lm_old = multipliers_[slave];
		lm_max += lm_old * lm_old;
		lm_old -= s2;
		lm_diff += lm_old * lm_old;
    
		// increase iterator
		++it;
	}
  
	if (lm_max < tol) {
		error = sqrt(lm_diff);
		return sqrt(lm_diff) < tol;
	}
  
	error = sqrt(lm_diff / lm_max);
	return sqrt(lm_diff / lm_max) < tol;
}

__END_AKANTU__

