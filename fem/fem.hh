/**
 * @file   fem.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jul 16 10:24:24 2010
 *
 * @brief  FEM class
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

#ifndef __AKANTU_FEM_HH__
#define __AKANTU_FEM_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "mesh.hh"
#include "element_class.hh"
#include "sparse_matrix.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__


class FEM : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FEM(Mesh & mesh, UInt spatial_dimension = 0,
      FEMID id = "fem", MemoryID memory_id = 0);

  virtual ~FEM();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// build the profile of the sparse matrix corresponding to the mesh
  void initSparseMatrixProfile(SparseMatrixType sparse_matrix_type = _unsymmetric);


  /// pre-compute all the shape functions, their derivatives and the jacobians
  virtual void initShapeFunctions(GhostType ghost_type = _not_ghost)=0;

  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */

  /// integrate f for all elements of type "type"
  virtual void integrate(const Vector<Real> & f,
		 Vector<Real> &intf,
		 UInt nb_degre_of_freedom,
   		 const ElementType & type,
   		 GhostType ghost_type = _not_ghost,
   		 const Vector<UInt> * filter_elements = NULL) const = 0;

  /// integrate a scalar value on all elements of type "type"
  virtual Real integrate(const Vector<Real> & f,
   		 const ElementType & type,
   		 GhostType ghost_type = _not_ghost,
   		 const Vector<UInt> * filter_elements = NULL) const = 0;

  /* ------------------------------------------------------------------------ */
  /* compatibility with old FEM fashion */
  /* ------------------------------------------------------------------------ */

  /// get the number of quadrature points
  virtual UInt getNbQuadraturePoints(const ElementType & type)=0;
  /// get the precomputed shapes
  const virtual Vector<Real> & getShapes(const ElementType & type)=0;
  /// get the precomputed ghost shapes
  const virtual Vector<Real> & getGhostShapes(const ElementType & type)=0;
  /// get the derivatives of shapes
  const virtual Vector<Real> & getShapesDerivatives(const ElementType & type)=0;
  /// get the ghost derivatives of shapes
  const virtual Vector<Real> & getGhostShapesDerivatives(const ElementType & type)=0;
  /// get quadrature points
  const virtual Vector<Real> & getQuadraturePoints(const ElementType & type)=0;

  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */

  virtual
  void gradientOnQuadraturePoints(const Vector<Real> &u,
				  Vector<Real> &nablauq,
				  const UInt nb_degre_of_freedom,
				  const ElementType & type,
				  GhostType ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL)=0;

  virtual
  void interpolateOnQuadraturePoints(const Vector<Real> &u,
				     Vector<Real> &uq,
				     UInt nb_degre_of_freedom,
				     const ElementType & type,
				     GhostType ghost_type = _not_ghost,
				     const Vector<UInt> * filter_elements = NULL) const =0;

  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */

  /// pre-compute normals on control points
  virtual void computeNormalsOnControlPoints(GhostType ghost_type = _not_ghost)=0;



  /// assemble vectors
  void assembleVector(const Vector<Real> & elementary_vect,
		      Vector<Real> & nodal_values,
		      UInt nb_degre_of_freedom,
		      const ElementType & type,
		      GhostType ghost_type = _not_ghost,
		      const Vector<UInt> * filter_elements = NULL,
		      Real scale_factor = 1) const;


  /// assemble matrix in the complete sparse matrix
  void assembleMatrix(const Vector<Real> & elementary_mat,
		      SparseMatrix & matrix,
		      const Vector<Int> & equation_number,
		      UInt nb_degre_of_freedom,
		      const ElementType & type,
		      GhostType ghost_type = _not_ghost,
		      const Vector<UInt> * filter_elements = NULL) const;


  /// assemble a field as a lumped matrix (ex. rho in lumped mass)
  virtual void assembleFieldLumped(__attribute__ ((unused)) const Vector<Real> & field_1,
				   __attribute__ ((unused)) Vector<Real> & lumped,
				   __attribute__ ((unused)) ElementType type,
				   __attribute__ ((unused)) __attribute__ ((unused)) GhostType ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  /// assemble a field as a matrix (ex. rho to mass matrix)
  virtual void assembleFieldMatrix(const Vector<Real> & field_1,
				   UInt nb_degree_of_freedom,
				   const Vector<Int> & equation_number,
				   SparseMatrix & matrix,
				   ElementType type,
				   GhostType ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }


  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// initialise the class
  void init();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(ElementDimension, element_dimension, UInt);

  /// get the mesh contained in the fem object
  inline Mesh & getMesh() const;

  /// get the in-radius of an element
  static inline Real getElementInradius(Real * coord, const ElementType & type);

  /// get the normals on quadrature points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(NormalsOnQuadPoints, normals_on_quad_points, const Vector<Real> &);


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// id of the fem object
  FEMID id;

  /// spatial dimension of the problem
  UInt element_dimension;

  /// the mesh on which all computation are made
  Mesh * mesh;

  /// normals at quadrature points
  ByElementTypeReal normals_on_quad_points;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "fem_inline_impl.cc"


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const FEM & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#include "fem_template.hh"

#endif /* __AKANTU_FEM_HH__ */
