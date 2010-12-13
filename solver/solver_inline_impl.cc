/**
 * @file   solver_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 17 16:14:58 2010
 *
 * @brief  implementation of inline function of the Solver class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
inline void Solver::assemble(const Vector<Real> & local_matrix, const Element & element) {
  AKANTU_DEBUG_ASSERT(matrix != NULL, "The sparse matrix must be first instantiated.");

  matrix->addToMatrix(local_matrix, element);
}
