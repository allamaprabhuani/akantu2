/**
 * @file   solver_mumps.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 17 17:28:56 2010
 *
 * @brief  Solver class implementation for the mumps solver
 *
 * @section LICENSE
 *
 * <insert license here>
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

  SolverMumps(const Mesh & mesh,
	      const SparseMatrixType & sparse_matrix_type,
	      UInt nb_degre_of_freedom,
	      const SolverID & id = "solver_mumps",
	      const MemoryID & memory_id = 0);

  virtual ~SolverMumps();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// build the profile and do the analysis part
  void initialize();

  /// factorize and solve the system
  void solve();

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;

private:

  inline Int & icntl(UInt i) {
    return mumps_data.icntl[i - 1];
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
