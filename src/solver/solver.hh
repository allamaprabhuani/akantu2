/**
 * @file   solver.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Oct  4 20:27:42 2010
 *
 * @brief  interface for solvers
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

#ifndef __AKANTU_SOLVER_HH__
#define __AKANTU_SOLVER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "aka_vector.hh"
#include "sparse_matrix.hh"
#include "mesh.hh"
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Solver : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Solver(SparseMatrix & matrix,
	 const ID & id = "solver",
	 const MemoryID & memory_id = 0);

  virtual ~Solver();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the solver
  virtual void initialize() = 0;

  /// solve
  virtual void solve(Vector<Real> & solution) = 0;
  virtual void solve() = 0;

  virtual void setRHS(Vector<Real> & rhs) = 0;

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(RHS, *rhs, Vector<Real> &);


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of thr Solver
  ID id;

  /// the matrix
  SparseMatrix * matrix;

  /// is sparse matrix allocated
  bool is_matrix_allocated;

  /// the right hand side
  Vector<Real> * rhs;

  /// mesh
  const Mesh * mesh;

  /// pointer to the communicator
  StaticCommunicator * communicator;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "solver_inline_impl.cc"
#endif

/// standard output stream operator
// __aka_inline__ std::ostream & operator <<(std::ostream & stream, const Solver & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__

#endif /* __AKANTU_SOLVER_HH__ */
