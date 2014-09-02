/**
 * @file   contact_manager_implicit.hh
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

#ifndef __AKANTU_IMPLICIT_CONTACT_MANAGER_HH__
#define __AKANTU_IMPLICIT_CONTACT_MANAGER_HH__

#include <unordered_set>

#include "mesh_graph.hh"
#include "model_manager.hh"
#include "dumper_iohelper.hh"
#include "aka_optimize.hh"
#include "integration_scheme_2nd_order.hh"


#define DEBUG_MANAGER 1


__BEGIN_AKANTU__

template <typename T>
T Heaviside(T v) {
	return v < 0 ? 0 : 1;
}

template <typename T>
T Macauley(T v) {
	return v < 0 ? 0 : v;
}

//! Enumerated type used for the ContactData overloaded operator[] that returns real values
enum Contact_parameter_type { Epsilon, Alpha, Uzawa_tol, Newton_tol, Uzawa_max_steps, Newton_max_steps};

//! Enumerated type used for the ContactData overloaded operator[] that returns boolean values
enum Contact_flag_type { Dump_iteration, Verbose, Automatic_penalty_parameter};


template <class vector_type, class point_type>
vector_type invert_map(const point_type& s, const point_type& p, const point_type& q) {
	vector_type b(2);

	b[0] = q[0] - p[0];
	b[1] = q[1] - p[1];

	vector_type a(2);

	a[0] = s[0] - p[0];
	a[1] = s[0] - p[1];

	vector_type r(1);

	Real s1 = transpose(a) * b;
	Real s2 = transpose(b) * b;
	r[0] = 2 * s1 / (s2) - 1;
	return r;
}

template <class vector_type>
void compute_shapes(const vector_type& natural_coords, vector_type& N) {
	/// natural coordinate
	Real c = natural_coords(0);
	/// shape functions
	N(0) = 0.5 * (1 - c);
	N(1) = 0.5 * (1 + c);
}

//! Function template specialization for inversion of a \f$ 3 \times 3 \f$ matrix.
template <class matrix_type>
std::pair <matrix_type, Real> invert(matrix_type& A) {
	// obtain determinant of the matrix
	Real det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
	    A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
	    A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

	// compute inverse
	matrix_type inv(3, 3, 1. / det);

	inv[0][0] *=  A[1][1] * A[2][2] - A[1][2] * A[2][1];
	inv[0][1] *=  A[0][2] * A[2][1] - A[0][1] * A[2][2];
	inv[0][2] *= -A[0][2] * A[1][1] + A[0][1] * A[1][2];
	inv[1][0] *=  A[1][2] * A[2][0] - A[1][0] * A[2][2];
	inv[1][1] *=  A[0][0] * A[2][2] - A[0][2] * A[2][0];
	inv[1][2] *=  A[0][2] * A[1][0] - A[0][0] * A[1][2];
	inv[2][0] *= -A[1][1] * A[2][0] + A[1][0] * A[2][1];
	inv[2][1] *=  A[0][1] * A[2][0] - A[0][0] * A[2][1];
	inv[2][2] *= -A[0][1] * A[1][0] + A[0][0] * A[1][1];

	return std::make_pair(inv, det);
}

template <class vector_type, class point_type>
vector_type invert_map(const point_type& s, const point_type& a, const point_type& b, const point_type& c) {
	typedef array::Array <2, Real> matrix_type;

	// matrix for inverse
	matrix_type A = {
		{ b[0] - a[0], c[0] - a[0], a[0] },
		{ b[1] - a[1], c[1] - a[1], a[1] },
		{ b[2] - a[2], c[2] - a[2], a[2] }
	};

	std::pair <matrix_type, Real> Ainv = invert(A);
	vector_type x = { s[0], s[1], s[2] };
	vector_type r1 = Ainv.first * x;

	return vector_type { r1[0], r1[1] }; // return only the first two components of r1
}

template <class vector_type>
void compute_shape_derivatives(const vector_type&, vector_type& DN) {
	DN(0) = -0.5;
	DN(1) =  0.5;
}

template <int dim, class model_type>
class SearchTraits;

template <class model_type>
struct SearchTraits <3, model_type> {
	typedef Point <3> point_type;
	typedef ModelElement <model_type> element_type;

	static bool check_projection(const point_type& p, model_type& model, UInt id) {
		element_type m(model, _triangle_3, id);

		return point_has_projection_to_triangle(p, m.template point <3>(0), m.template point <3>(1), m.template point <3>(2));
	}
};


template <class model_type>
struct SearchTraits <2, model_type> {
	typedef Point <2> point_type;
	typedef ModelElement <model_type> element_type;

	static bool check_projection(const point_type& p, model_type& model, UInt id) {
		element_type m(model, _segment_2, id);

		return has_projection(p, m.template point <2>(0), m.template point <2>(1));
	}
};



template <int dim, class model_type>
struct ContactData {
	typedef Point <dim> point_type;
	typedef array::Array <1, Real> vector_type;
	typedef array::Array <2, Real> matrix_type;
	typedef ModelElement <model_type> element_type;

	typedef std::map <Contact_parameter_type, Real> options_map;
	typedef std::set <Contact_flag_type> flag_map;

	typedef std::map <UInt, element_type> slave_master_map;
	typedef std::map <UInt, Real> real_map;
	typedef typename real_map::iterator real_iterator;
	typedef std::map <UInt, Real> gap_map;

	struct SearchBase {
		virtual int search(const Real *) {
			static bool msg = false;
			if (!msg) {
				std::cout << " - Warning: calling default base searcher, type any key to continue." << std::endl;
				std::cin.ignore();
				msg = true;
			}
			return -1;
		}

		virtual ~SearchBase()
		{
		}
	};

	struct MasterAssignator : public SearchBase {
		std::string surface_;
		model_type &model_;

		MasterAssignator(const std::string & s, model_type & model) : surface_(s), model_(model)
		{
		}

		virtual int search(const Real *ptr) {
			point_type p(ptr);

			ElementGroup &rs = model_.getMesh().getElementGroup(surface_);

			for (ElementGroup::type_iterator tit = rs.firstType(); tit != rs.lastType(); ++tit)
				for (ElementGroup::const_element_iterator it = rs.element_begin(*tit);
				     it != rs.element_end(*tit); ++it) {
					if (SearchTraits <dim, model_type>::check_projection(p, model_, *it))
						return *it;
				}
			return -1;
		}
	};

	typedef SearchBase search_type;


	slave_master_map sm_;
	real_map multipliers_, areas_, penalty_, gaps_;
	model_type &model_;
	Array <Real> multiplier_dumper_, pressure_dumper_;
	options_map options_;
	flag_map flags_;
	size_t uiter_, niter_;
	SearchBase *searcher_;

	ContactData(int argc, char *argv[], model_type & m) : model_(m), multiplier_dumper_(m.getMesh().getNbNodes(), 3), pressure_dumper_(m.getMesh().getNbNodes(), 3), options_(), uiter_(), niter_(), searcher_(new SearchBase())
	{
		// register dumpers
   	        model_.getMesh().addDumpFieldExternal("multipliers", multiplier_dumper_);
		model_.getMesh().addDumpFieldExternal("pressure", pressure_dumper_);

		set_default_options();

		cout << "\nValid parameters for contact options:" << endl;
		cout << "  -e      [auto] Penalty parameter for Augmented-Lagrangian formulation" << endl;
		cout << "  -alpha  [1] Multiplier for values of the penalty parameter" << endl;
		cout << "  -utol   [1e-4] Tolerance used for multipliers in the Uzawa method" << endl;
		cout << "  -ntol   [1e-4] Tolerance used in the Newton-Raphson inner convergence loop" << endl;
		cout << "  -usteps [100]  Maximum number of steps allowed in the Uzawa loop" << endl;
		cout << "  -nsteps [100]  Maximum number of steps allowed in the Newton-Raphson loop" << endl;
		cout << "\nValid flags:" << endl;
		cout << "  -dump          Dumping within Newton iterations" << endl;
		cout << "  -v             Verbose flag\n" << endl;

		flags_.insert(Automatic_penalty_parameter);

		options_[Alpha] = 1.;

		// loop over command line parameters
		for (int i = 0; i < argc; ++i) {
			if (strcmp(argv[i], "-e") == 0) {
				Real epsilon = atof(argv[++i]);
				options_[Epsilon] = epsilon;
				cout << "-e = " << epsilon << endl;
				cout << "Removing automatic penalty parameter computation" << endl;
				flags_.erase(Automatic_penalty_parameter);
			}
			else if (strcmp(argv[i], "-alpha") == 0) {
				Real alpha = atof(argv[++i]);
				assert(alpha > 0);
				options_[Alpha] = alpha;
				cout << "-alpha = " << alpha << endl;
			}
			else if (strcmp(argv[i], "-utol") == 0) {
				Real utol = atof(argv[++i]);
				assert(utol > 0);
				options_[Uzawa_tol] = utol;
				cout << "-utol = " << utol << endl;
			}
			else if (strcmp(argv[i], "-ntol") == 0) {
				Real ntol = atof(argv[++i]);
				assert(ntol > 0);
				options_[Newton_tol] = ntol;
				cout << "-ntol = " << ntol << endl;
			}
			else if (strcmp(argv[i], "-usteps") == 0) {
				Real usteps = atof(argv[++i]);
				assert(usteps > 0);
				options_[Uzawa_max_steps] = usteps;
				cout << "-usteps = " << usteps << endl;
			}
			else if (strcmp(argv[i], "-nsteps") == 0) {
				Real nsteps = atof(argv[++i]);
				assert(nsteps > 0);
				options_[Newton_max_steps] = nsteps;
				cout << "-nsteps = " << nsteps << endl;
			}
			else if (strcmp(argv[i], "-dump") == 0) {
				cout << "-dump flag given" << endl;
				flags_.insert(Dump_iteration);
			}
			else if (strcmp(argv[i], "-eauto") == 0) {
				cout << "-eauto flag given" << endl;
				flags_.insert(Automatic_penalty_parameter);
			}
			else if (strcmp(argv[i], "-v") == 0) {
				cout << "-v flag given" << endl;
				flags_.insert(Verbose);
			}
		}
	}


	struct PostAssemblyEmptyFunctor {
		void operator()() const {
		}
	};

	template <class PostAssemblyFunctor>
	void solveContactCommon(SearchBase *sf = new SearchBase(), const PostAssemblyFunctor& paf = PostAssemblyEmptyFunctor()) {
		ContactData &cd = *this;

		model_.implicitPred();
		model_.updateResidual();

		AKANTU_DEBUG_ASSERT(model_.stiffness_matrix != NULL,
		                    "You should first initialize the implicit solver and assemble the stiffness matrix");


		// implementation of the Uzawa method for solving contact
		bool uzawa_converged = false;
		static UInt step = 0;
		UInt k = 0;
		UInt ntotal = 0;

		std::list <int> nt;


		std::ofstream ofs;
		ofs.open("iterations.out", std::ofstream::out | std::ofstream::app);

		// initialize Lagrange multipliers
		// NOTE: It doesn't make any difference to start from the previous
		// converged solution of Lagrange multipliers
		real_map lambda_new;

		cout << std::boolalpha;

		if (cd[Verbose])
			cout << "- Start Uzawa:" << endl;

		do {
			Real uerror = 0.;

			bool converged = false;
			UInt j = 0;

			cd.uiter_ = k;

			do {
				Real nerror = 0.;

				cd.niter_ = j;

				// assemble material matrix
				model_.assembleStiffnessMatrix();

				// call post-assembly functor
				paf();

				// compute gaps
				uzawa_converged = cd.computeTangentAndResidual(lambda_new, sf, uerror);

				// solve
				model_.template solve <IntegrationScheme2ndOrder::_displacement_corrector> (*model_.increment);

				model_.implicitCorr();
				model_.updateResidual();

				converged = model_.template testConvergence <_scc_increment> (cd[Newton_tol], nerror);

				if (cd[Dump_iteration])
					cd.dump();
				if (cd[Verbose])
					cout << "    Newton: " << j << ", " << nerror << " < " << cd[Newton_tol] << " = " << (nerror < cd[Newton_tol]) << endl;

				++j;
				AKANTU_DEBUG_INFO("[" << _scc_increment << "] Convergence iteration "
				                      << std::setw(std::log10(cd[Newton_max_steps])) << j
				                      << ": error " << nerror << (converged ? " < " : " > ") << cd[Newton_tol] << std::endl);
			}
			while (!converged && j < cd[Newton_max_steps]);


			if (cd[Verbose])
				cout << "  Uzawa: " << k << ", " << uerror << " < " << cd[Uzawa_tol] << " = " << (uerror < cd[Uzawa_tol]) << endl;

			if (j == cd[Newton_max_steps]) {
				cout << "*** ERROR *** Newton-Raphson loop did not converge within max number of iterations: " << cd[Newton_max_steps] << endl;
				exit(1);
			}

			nt.push_back(j);
			ntotal += j;

			// increment uzawa loop counter
			++k;

			AKANTU_DEBUG_INFO("[" << _scc_increment << "] Uzawa convergence iteration "
			                      << std::setw(std::log10(cd[Newton_max_steps])) << k
			                      << std::endl);

			// update lagrange multipliers
			cd.multipliers_ = lambda_new;
		}
		while (!uzawa_converged && k < cd[Uzawa_max_steps]);

		if (k == cd[Uzawa_max_steps]) {
			cout << "*** ERROR *** Uzawa loop did not converge within max number of iterations: " << cd[Uzawa_max_steps] << endl;
			exit(1);
		}


		cout << "Summary: Uzawa [" << k << "]: Newton [" << ntotal << "]:";
		for (int n : nt)
			cout << " " << n;
		cout << endl;

		ofs << std::setw(10) << ++step << std::setw(10) << k << std::setw(10) << ntotal << endl;
		ofs.close();
	}

	void searchSurface(const std::string& s) {
		delete searcher_;
		searcher_ = new MasterAssignator(s, model_);
	}

	void solveContactStep() {
		solveContactCommon <PostAssemblyEmptyFunctor>(searcher_);

		dump();
	}

	template <class PostAssemblyFunctor>
	void solveContactStep(const PostAssemblyFunctor& fn) {
		solveContactCommon <PostAssemblyFunctor>(searcher_, fn);

		dump();
		fn.dump();
	}

//	void solveContactStep() {
//		solveContactCommon <PostAssemblyEmptyFunctor>();
//
//		dump();
//	}
//
//	void solveContactStep(search_type *search) {
//		solveContactCommon <PostAssemblyEmptyFunctor>(search);
//		dump();
//	}
//
//	template <class PostAssemblyFunctor>
//	void solveContactStep(search_type *search, const PostAssemblyFunctor& fn) {
//		solveContactCommon <PostAssemblyFunctor>(search, fn);
//
//		dump();
//		fn.dump();
//	}

	void prepareDump() {
		multiplier_dumper_.clear();
		pressure_dumper_.clear();

		for (auto v : multipliers_) {
			element_type &el = sm_[v.first];
			auto n = el.normal();

			Real lambda = v.second;

			for (int i = 0; i < n.size(); ++i)
				multiplier_dumper_(v.first, i) = lambda * n[i];

			// dump pressures only if area is associated with node
			auto it = areas_.find(v.first);
			if (it != areas_.end())
				for (int i = 0; i < n.size(); ++i) {
					Real a = it->second;
					assert(a != 0.);
					pressure_dumper_(v.first, i) = lambda * n[i] / a;
				}
			else
				cout << "*** WARNING *** Zero area for slave node " << v.first << endl;
		}
	}

	//! Function that dumps paraview files
	void dump() {
		prepareDump();
		model_.dump();
	}

	//! Overloaded operator[] (const) that returns real values
	Real operator[] (Contact_parameter_type p)const {
		auto it = options_.find(p);
		assert(it != options_.end());
		return it->second;
	}

	//! Overloaded operator[] (non-const) that returns real values
	Real& operator[] (Contact_parameter_type p)
	{return options_[p]; }

	//! Overloaded operator[] (const) that returns flags
	bool operator[] (Contact_flag_type f)const {
		auto it = flags_.find(f);
		return it != flags_.end();
	}

	//! Get contact force (sum of Lagrange multipliers)
	Real getForce() {
		Real f = 0.;
		for (auto v : multipliers_)
			f += v.second;
		return f;
	}

	//! Get maximum value of contact pressure
	Real getMaxPressure() {
		Real p = 0.;
		for (auto v : multipliers_) {
			auto slave = v.first;
			auto it = areas_.find(slave);
			if (it != areas_.end())
				p = std::max(p, v.second / it->second); // max(p, lambda/area)
		}
		return p;
	}

	//! Add slave
	void addSlave(UInt s) {
		sm_[s] = element_type();
		multipliers_[s] = Real();
	}

	//! Add slave-master pair
	void addPair(UInt s, element_type el) {
		sm_[s] = el;
		multipliers_[s] = Real();
	}

	//! Add area to a slave node
	void addArea(UInt n, Real a) {
		if (a != 0.) areas_[n] = a;
	}

	//! Compute contact contributions to stiffness matrix and force vector
	bool computeTangentAndResidual(real_map&, search_type *, Real&);

	//! Provide standard output of contact object
	friend std::ostream& operator << (std::ostream & os, const ContactData &cd) {
		os << "\n*** INFO *** ContactData info:" << endl;
		os << "Contact parameters:" << endl;
		if (cd[Automatic_penalty_parameter])
			cout << "  -e      = auto" << endl;
		else
			cout << "  -e      = " << cd[Epsilon] << endl;

		cout << "  -alpha  = " << cd[Alpha] << endl;
		cout << "  -utol   = " << cd[Uzawa_tol] << endl;
		cout << "  -ntol   = " << cd[Newton_tol] << endl;

		cout << "  -usteps = " << cd[Uzawa_max_steps] << endl;
		cout << "  -nsteps = " << cd[Newton_max_steps] << endl;

		cout << "\nContact flags:" << endl;
		cout << "  -dump   = " << cd[Dump_iteration] << endl;
		cout << "  -v      = " << cd[Verbose] << endl;

		if (cd[Verbose]) {
			// loop over pairs
			cout << "\nSlave nodes: ";
			for (auto it = cd.sm_.begin(); it != cd.sm_.end(); ++it)
				os << it->first << " ";
			os << endl;

			// loop over pairs
			cout << "\nSlave master pairs" << endl;
			for (auto it = cd.sm_.begin(); it != cd.sm_.end(); ++it) {
				auto slave = it->first;
				auto master = it->second;
				os << "  slave: " << slave << ", Master: ";
				if (master == element_type())
					os << "none" << endl;
				else
					os << master << endl;
			}
		}
		return os;
		}

	    private:
	    //! Compute penalty parameter values for each slave node automatically
	    void getPenaltyValues() {
		cout << "*** INFO *** Obtaining penalty parameters automatically. ";
		const SparseMatrix &Kconst = model_.getStiffnessMatrix();

		Real ave = 0.;
		size_t k = 0;

		// loop over pairs
		for (auto it = sm_.begin(); it != sm_.end(); ++it) {
			auto slave = it->first;
			auto master = it->second;

			if (master != element_type()) {
				std::vector <UInt> conn(master.numNodes() + 1); // 1 slave (not hardcoded)

				conn[0] = slave;
				for (UInt i = 0; i < master.numNodes(); ++i)
					conn[1 + i] = master.node(i);

				// compute normal
				vector_type nu = master.normal();

				// carry out stiffness multiplication with the normal
				// the product Kij*nj would give the force for a unit displacement
				// (i.e., the stiffness needed to move the node by 1)
				matrix_type r(Kconst.getSize(), master.numNodes() + 1);

				// loop over stifness matrix dimension
				for (int i = 0; i < Kconst.getSize(); ++i)
					// loop over problem dimensions
					for (int j = 0; j < dim; ++j)
						// loop over nodes considered
						for (int k = 0; k < master.numNodes() + 1; ++k)
							r(i, k) += Kconst(i, conn[k] + j) * nu(j);

				// get results (norm of each column in r)
				vector_type rsum(master.numNodes() + 1);

				for (int i = 0; i < rsum.size(); ++i)
					for (int j = 0; j < r.rows(); ++j)
						rsum(i) += r(j, i) * r(j, i);

				// get average value as the penalty parameter
				Real epsilon = 0.;
				for (int i = 0; i < rsum.size(); ++i)
					epsilon += sqrt(rsum(i));

				epsilon /= master.numNodes() + 1;
				penalty_[slave] = epsilon;

				ave += penalty_[slave];
				++k;
			}
			// dummy master
			else {
				// carry out stiffness multiplication with the normal
				// the product Kij*nj would give the force for a unit displacement
				// (i.e., the stiffness needed to move the node by 1)
				vector_type r(Kconst.getSize());

				// loop over stifness matrix dimension
				for (int i = 0; i < Kconst.getSize(); ++i)
					// loop over problem dimensions
					for (int j = 0; j < dim; ++j)
						// loop over nodes considered
						r(i) += Kconst(i, slave + j) * 1. / dim;

				// get results (norm of each column in r)
				Real epsilon = 0;
				for (int i = 0; i < r.size(); ++i)
					epsilon += r(i) * r(i);

				epsilon = sqrt(epsilon);
				penalty_[slave] = epsilon;

				ave += penalty_[slave];
				++k;
			}
		}
		cout << "Average value: " << (*this)[Alpha] * ave / k << endl;
		}

	    // set default options
	    void set_default_options() {
		options_[Epsilon] = 1000;
		options_[Uzawa_tol] = 1e-4;
		options_[Newton_tol] = 1e-4;
		options_[Uzawa_max_steps] = 100;
		options_[Newton_max_steps] = 100;
		}
};


__END_AKANTU__

#endif /* __AKANTU_IMPLICIT_CONTACT_MANAGER_HH__ */
