/**
 * @file   sparse_matrix.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Oct  4 20:33:00 2010
 *
 * @brief  sparse matrix storage class (distributed assembled matrix)
 * This is a COO format (Coordinate List)
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

  // SparseMatrix(const Mesh & mesh,
  // 	       const SparseMatrixType & sparse_matrix_type,
  // 	       UInt nb_degre_of_freedom,
  // 	       const SparseMatrixID & id = "sparse_matrix",
  // 	       const MemoryID & memory_id = 0);

  SparseMatrix(UInt size,
	       const SparseMatrixType & sparse_matrix_type,
	       UInt nb_degre_of_freedom,
	       const SparseMatrixID & id = "sparse_matrix",
	       const MemoryID & memory_id = 0);

  SparseMatrix(const SparseMatrix & matrix);

  virtual ~SparseMatrix();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// remove the existing profile
  inline void clearProfile();

  /// add a non-zero element
  UInt addToProfile(UInt i, UInt j);

  /// set the matrix to 0
  inline void clear();

  /// assemble a local matrix in the sparse one
  inline void addToMatrix(UInt i, UInt j, Real value);

  /// fill the profil of the matrix
  void buildProfile(const Mesh & mesh, const Vector<Int> & equation_number);

  /// modify the matrix to "remove" the blocked dof
  void applyBoundary(const Vector<bool> & boundary,
		     const unordered_map<UInt, UInt>::type & local_eq_num_to_global);

  /// modify the matrix to "remove" the blocked dof
  void removeBoundary(const Vector<bool> & boundary);

  /// restore the profile that was before removing the boundaries
  void restoreProfile();

  /// save the profil in a file using the MatrixMarket file format
  void saveProfile(const std::string & filename);

  /// save the matrix in a file using the MatrixMarket file format
  void saveMatrix(const std::string & filename);

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  inline UInt key(UInt i, UInt j) const {
    if(sparse_matrix_type == _symmetric && (i > j))
      return j * size + i;

    return i * size + j;
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the values at potition i, j
  inline Real operator()(UInt i, UInt j) const;
  inline Real & operator()(UInt i, UInt j);

  AKANTU_GET_MACRO(IRN, irn, const Vector<Int> &);

  AKANTU_GET_MACRO(JCN, jcn, const Vector<Int> &);

  AKANTU_GET_MACRO(A, a, const Vector<Real> &);

  AKANTU_GET_MACRO(NbNonZero, nb_non_zero, UInt);

  AKANTU_GET_MACRO(Size, size, UInt);

  AKANTU_GET_MACRO(NbDegreOfFreedom, nb_degre_of_freedom, UInt);

  AKANTU_GET_MACRO(SparseMatrixType, sparse_matrix_type, const SparseMatrixType &);

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
  //  const Mesh * mesh;

  /// Size of the matrix
  UInt size;

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


  /// saved row indexes
  Vector<Int> * irn_save;

  /// saved column indexes
  Vector<Int> * jcn_save;

  /// saved size
  UInt size_save;

  /// information to know where to assemble an element in a global sparse matrix
  //  ByElementTypeUInt element_to_sparse_profile;

  /* map for  (i,j) ->  k correspondence \warning std::map are slow
   *  \todo improve  with hash_map (non standard in stl) or unordered_map (boost or C++0x)
   */
  typedef unordered_map<UInt, UInt>::type coordinate_list_map;
  coordinate_list_map irn_jcn_k;

  std::map<std::pair<UInt, UInt>, UInt> * irn_jcn_to_k;
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

Vector<Real> & operator*=(Vector<Real> & vect, const SparseMatrix & mat);


__END_AKANTU__

#endif /* __AKANTU_SPARSE_MATRIX_HH__ */
