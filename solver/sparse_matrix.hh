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
#ifndef __INTEL_COMPILER
namespace std {
  namespace tr1 {
    template<typename a, typename b>
    struct hash< std::pair<a, b> > {
    private:
      const hash<a> ah;
      const hash<b> bh;
    public:
      hash() : ah(), bh() {}
      size_t operator()(const std::pair<a, b> &p) const {
	size_t seed = ah(p.first);
	return bh(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
      }
    };
  }
} // namespaces
#endif

__BEGIN_AKANTU__

class DOFSynchronizer;

class SparseMatrix : private Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrix(UInt size,
	       const SparseMatrixType & sparse_matrix_type,
	       UInt nb_degre_of_freedom,
	       const SparseMatrixID & id = "sparse_matrix",
	       const MemoryID & memory_id = 0);

  SparseMatrix(const SparseMatrix & matrix,
	       const SparseMatrixID & id = "sparse_matrix",
	       const MemoryID & memory_id = 0);

  virtual ~SparseMatrix();

  typedef std::pair<UInt, UInt> KeyCOO;
  typedef unordered_map<KeyCOO, UInt>::type coordinate_list_map;

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
  void buildProfile(const Mesh & mesh, const DOFSynchronizer & dof_synchronizer);

  /// modify the matrix to "remove" the blocked dof
  void applyBoundary(const Vector<bool> & boundary);

  /// modify the matrix to "remove" the blocked dof
  void removeBoundary(const Vector<bool> & boundary);

  /// restore the profile that was before removing the boundaries
  void restoreProfile();

  /// save the profil in a file using the MatrixMarket file format
  void saveProfile(const std::string & filename);

  /// save the matrix in a file using the MatrixMarket file format
  void saveMatrix(const std::string & filename);

  /// copy assuming the profil are the same
  void copyContent(const SparseMatrix & matrix);

  /// add matrix assuming the profil are the same
  void add(const SparseMatrix & matrix, Real alpha);

  /// diagonal lumping
  void lump(Vector<Real> & lumped);

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  inline KeyCOO key(UInt i, UInt j) const {
    if(sparse_matrix_type == _symmetric && (i > j))
      return std::make_pair(j, i);

    return std::make_pair(i, j);
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

  const DOFSynchronizer & getDOFSynchronizer() const {
    AKANTU_DEBUG_ASSERT(dof_synchronizer != NULL,
			"DOFSynchronizer not initialized in the SparseMatrix!");
    return *dof_synchronizer;
  }

private:
  AKANTU_GET_MACRO(DOFSynchronizerPointer, dof_synchronizer, DOFSynchronizer *);

  friend Vector<Real> & operator*=(Vector<Real> & vect, const SparseMatrix & mat);
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
  coordinate_list_map irn_jcn_k;

  DOFSynchronizer * dof_synchronizer;
  //  std::map<std::pair<UInt, UInt>, UInt> * irn_jcn_to_k;
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
