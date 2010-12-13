/**
 * @file   sparse_matrix.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Oct  4 20:33:00 2010
 *
 * @brief  sparse matrix storage class (distributed assembled matrix)
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SPARSE_MATRIX_HH__
#define __AKANTU_SPARSE_MATRIX_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class SparseMatrix : private Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SparseMatrix(const Mesh & mesh,
	       const SparseMatrixType & sparse_matrix_type,
	       UInt nb_degre_of_freedom,
	       const SparseMatrixID & id = "sparse_matrix",
	       const MemoryID & memory_id = 0);
  //  virtual ~SparseMatrix();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// add a non-zero element
  UInt addToProfile(UInt i, UInt j);

  // /// add a non-zero element with local number of nodes
  // UInt addToProfileLocal(UInt i, UInt j);

  /// assemble a local matrix in the sparse one
  inline void addToMatrix(const Vector<Real> & local_matrix,
			  const Element & element);

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

  /// fill the profil of the matrix
  void buildProfile();

  /// save the profil in a file using the MatrixMarket file format
  void saveProfile(const std::string & filename);

  /// save the matrix in a file using the MatrixMarket file format
  void saveMatrix(const std::string & filename);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(IRN, irn, const Vector<Int> &);

  AKANTU_GET_MACRO(JCN, jcn, const Vector<Int> &);

  AKANTU_GET_MACRO(A, a, const Vector<Real> &);

  AKANTU_GET_MACRO(NbNonZero, nb_non_zero, UInt);

  AKANTU_GET_MACRO(NbDegreOfFreedom, nb_degre_of_freedom, UInt);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// id  of the SparseMatrix
  SparseMatrixID id;

  /// sparce matrix type
  SparseMatrixType sparse_matrix_type;

  /// number of degre of freedom
  UInt nb_degre_of_freedom;

  /// Mesh corresponding to the profile
  const Mesh & mesh;

  /// number of processors
  UInt nb_proc;

  /// number of non zero element
  UInt nb_non_zero;

  /// row indexes
  Vector<Int> irn;

  /// column indexes
  Vector<Int> jcn;

  /// values : A[k] = Matrix[irn[k]][jcn[k]]
  Vector<Real> a;

  /// information to know where to assemble an element in a global sparse matrix
  ByElementTypeUInt element_to_sparse_profile;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "sparse_matrix_inline_impl.cc"

// /// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const SparseMatrix & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_SPARSE_MATRIX_HH__ */
