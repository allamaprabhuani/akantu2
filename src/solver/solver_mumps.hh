/**
 * @file   solver_mumps.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 17 17:28:56 2010
 *
 * @brief  Solver class implementation for the mumps solver
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


#ifndef __AKANTU_SOLVER_MUMPS_HH__
#define __AKANTU_SOLVER_MUMPS_HH__
#include <dmumps_c.h>

#include "solver.hh"

__BEGIN_AKANTU__

class SolverMumps : public Solver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SolverMumps(SparseMatrix & sparse_matrix,
	      const ID & id = "solver_mumps",
	      const MemoryID & memory_id = 0);

  virtual ~SolverMumps();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// build the profile and do the analysis part
  void initialize();

  /// factorize and solve the system
  void solve(Vector<Real> & solution);
  void solve();

  virtual void setRHS(Vector<Real> & rhs);

  // void initNodesLocation(const Mesh & mesh, UInt nb_degree_of_freedom);

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;

private:

  inline Int & icntl(UInt i) {
    return mumps_data.icntl[i - 1];
  }

  inline Int & info(UInt i) {
    return mumps_data.info[i - 1];
  }


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// mumps data
  DMUMPS_STRUC_C mumps_data;

  /// Vector containing the proc number for each node
  Vector<UInt> ** rhs_position;

  /// solution repartition over processors
  Vector<UInt> ** solution_position;

  /// nb node per proc
  UInt * nb_nodes_per_proc;

  /// nb node per proc
  UInt * nb_nodes_per_proc_rhs;

  /// nb local nodes
  UInt nb_local_nodes;

  /// nodes types
  Int * nodes_type;

  /* ------------------------------------------------------------------------ */
  /* Local types                                                              */
  /* ------------------------------------------------------------------------ */
private:

  bool rhs_is_local;

  enum SolverMumpsJob {
    _smj_initialize = -1,
    _smj_analyze = 1,
    _smj_factorize = 2,
    _smj_solve = 3,
    _smj_analyze_factorize = 4,
    _smj_factorize_solve = 5,
    _smj_complete = 6, // analyze, factorize, solve
    _smj_destroy = -2
  };
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "solver_mumps_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const SolverMumps & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_SOLVER_MUMPS_HH__ */
