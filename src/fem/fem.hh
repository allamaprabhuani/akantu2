/**
 * @file   fem.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Jul 20 23:40:43 2010
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
namespace akantu {
  class Integrator;
  class ShapeFunctions;
}

__BEGIN_AKANTU__

class QuadraturePoint : public Element {
public:
  typedef Vector<Real> position_type;
public:
  QuadraturePoint(ElementType type = _not_defined, UInt element = 0,
		  UInt num_point = 0, GhostType ghost_type = _not_ghost) :
    Element(type, element, ghost_type), num_point(num_point), global_num(0),
    position((Real *)NULL, 0) { };

  QuadraturePoint(UInt element, UInt num_point,
		  UInt global_num,
		  const position_type & position,
		  ElementType type,
		  GhostType ghost_type = _not_ghost) :
    Element(type, element, ghost_type), num_point(num_point), global_num(global_num),
    position((Real *)NULL, 0) { this->position.shallowCopy(position); };

  QuadraturePoint(const QuadraturePoint & quad) :
    Element(quad), num_point(quad.num_point), global_num(quad.global_num), position((Real *) NULL, 0) {
    position.shallowCopy(quad.position);
  };

  inline QuadraturePoint & operator=(const QuadraturePoint & q) {
    if(this != &q) {
      element    = q.element;
      type       = q.type;
      ghost_type = q.ghost_type;
      num_point  = q.num_point;
      global_num = q.global_num;
      position.shallowCopy(q.position);
    }

    return *this;
  }

  AKANTU_GET_MACRO(Position, position, const position_type &);

  void setPosition(const position_type & position) {
    this->position.shallowCopy(position);
  }

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
    stream << space << "QuadraturePoint [";
    Element::printself(stream, 0);
    stream << ", " << num_point << "]";
  }

public:
  UInt num_point;
  UInt global_num;
private:
  position_type position;
};

/**
 * The  generic  FEM class  derived  in  a  FEMTemplate class  containing  the
 * shape functions and the integration method
 */
class FEM : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FEM(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
      ID id = "fem", MemoryID memory_id = 0);

  virtual ~FEM();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// build the profile of the sparse matrix corresponding to the mesh
  void initSparseMatrixProfile(SparseMatrixType sparse_matrix_type = _unsymmetric);

  /// pre-compute all the shape functions, their derivatives and the jacobians
  virtual void initShapeFunctions(const GhostType & ghost_type = _not_ghost) = 0;

  /// extract the nodal values and store them per element
  template<typename T>
  static void extractNodalToElementField(const Mesh & mesh,
					 const Array<T> & nodal_f,
					 Array<T> & elemental_f,
					 const ElementType & type,
					 const GhostType & ghost_type = _not_ghost,
					 const Array<UInt> & filter_elements = empty_filter);

  /// filter a field
  template<typename T>
  static void filterElementalData(const Mesh & mesh,
                                  const Array<T> & quad_f,
                                  Array<T> & filtered_f,
                                  const ElementType & type,
                                  const GhostType & ghost_type = _not_ghost,
                                  const Array<UInt> & filter_elements = empty_filter);



  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */
  /// integrate f for all elements of type "type"
  virtual void integrate(const Array<Real> & f,
		 Array<Real> &intf,
		 UInt nb_degree_of_freedom,
   		 const ElementType & type,
   		 const GhostType & ghost_type = _not_ghost,
   		 const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// integrate a scalar value on all elements of type "type"
  virtual Real integrate(const Array<Real> & f,
   		 const ElementType & type,
   		 const GhostType & ghost_type = _not_ghost,
   		 const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// integrate f for all quadrature points of type "type"
  virtual void integrateOnQuadraturePoints(const Array<Real> & f,
					   Array<Real> &intf,
					   UInt nb_degree_of_freedom,
					   const ElementType & type,
					   const GhostType & ghost_type = _not_ghost,
					   const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// integrate one element scalar value on all elements of type "type"
  virtual Real integrate(const Vector<Real> & f,
			 const ElementType & type,
			 UInt index, const GhostType & ghost_type = _not_ghost) const = 0;


  /* ------------------------------------------------------------------------ */
  /* compatibility with old FEM fashion */
  /* ------------------------------------------------------------------------ */
  /// get the number of quadrature points
  virtual UInt getNbQuadraturePoints(const ElementType & type,
				     const GhostType & ghost_type = _not_ghost) const = 0;
  /// get the precomputed shapes
  const virtual Array<Real> & getShapes(const ElementType & type,
					const GhostType & ghost_type = _not_ghost) const = 0;

  /// get the derivatives of shapes
  const virtual Array<Real> & getShapesDerivatives(const ElementType & type,
						   const GhostType & ghost_type = _not_ghost,
						   UInt id = 0) const = 0;

  /// get quadrature points
  const virtual Matrix<Real> & getQuadraturePoints(const ElementType & type,
						   const GhostType & ghost_type = _not_ghost) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */
  virtual
  void gradientOnQuadraturePoints(const Array<Real> &u,
				  Array<Real> &nablauq,
				  const UInt nb_degree_of_freedom,
				  const ElementType & type,
				  const GhostType & ghost_type = _not_ghost,
                                  const Array<UInt> & filter_elements = empty_filter) const = 0;

  virtual
  void interpolateOnQuadraturePoints(const Array<Real> &u,
				     Array<Real> &uq,
				     UInt nb_degree_of_freedom,
				     const ElementType & type,
				     const GhostType & ghost_type = _not_ghost,
                                     const Array<UInt> & filter_elements = empty_filter) const = 0;

  virtual
  void interpolateOnQuadraturePoints(const Array<Real> & u,
				     ByElementTypeReal & uq,
                                     const ByElementTypeUInt * filter_elements = NULL) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */

  /// pre-compute normals on control points
  virtual void computeNormalsOnControlPoints(const GhostType & ghost_type = _not_ghost) = 0;

  /// pre-compute normals on control points
  virtual void computeNormalsOnControlPoints(__attribute__((unused)) const Array<Real> & field,
					     __attribute__((unused)) const GhostType & ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// pre-compute normals on control points
  virtual void computeNormalsOnControlPoints(__attribute__((unused)) const Array<Real> & field,
					     __attribute__((unused)) Array<Real> & normal,
					     __attribute__((unused)) const ElementType & type,
					     __attribute__((unused)) const GhostType & ghost_type = _not_ghost) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }




  /// assemble vectors
  void assembleArray(const Array<Real> & elementary_vect,
		      Array<Real> & nodal_values,
		      const Array<Int> & equation_number,
		      UInt nb_degree_of_freedom,
		      const ElementType & type,
		      const GhostType & ghost_type = _not_ghost,
		      const Array<UInt> & filter_elements = empty_filter,
		      Real scale_factor = 1) const;

  /// assemble matrix in the complete sparse matrix
  void assembleMatrix(const Array<Real> & elementary_mat,
		      SparseMatrix & matrix,
		      UInt nb_degree_of_freedom,
		      const ElementType & type,
		      const GhostType & ghost_type = _not_ghost,
		      const Array<UInt> & filter_elements = empty_filter) const;


  /// assemble a field as a lumped matrix (ex. rho in lumped mass)
  virtual void assembleFieldLumped(__attribute__ ((unused)) const Array<Real> & field_1,
				   __attribute__ ((unused)) UInt nb_degree_of_freedom,
  				   __attribute__ ((unused)) Array<Real> & lumped,
  				   __attribute__ ((unused)) const Array<Int> & equation_number,
  				   __attribute__ ((unused)) ElementType type,
  				   __attribute__ ((unused)) const GhostType & ghost_type) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };


  /// assemble a field as a matrix (ex. rho to mass matrix)
  virtual void assembleFieldMatrix(__attribute__ ((unused)) const Array<Real> & field_1,
				   __attribute__ ((unused)) UInt nb_degree_of_freedom,
				   __attribute__ ((unused)) SparseMatrix & matrix,
				   __attribute__ ((unused)) ElementType type,
				   __attribute__ ((unused)) const GhostType & ghost_type) const {
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
  static inline Real getElementInradius(const Matrix<Real> & coord, const ElementType & type);

  /// get the normals on quadrature points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(NormalsOnQuadPoints, normals_on_quad_points, Real);

  /// get cohesive element type for a given facet type
  static inline ElementType getCohesiveElementType(const ElementType & type_facet);

  virtual const ShapeFunctions & getShapeFunctionsInterface() const = 0;
  virtual const Integrator & getIntegratorInterface() const = 0;

  static inline InterpolationType getInterpolationType(const ElementType & el_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// spatial dimension of the problem
  UInt element_dimension;

  /// the mesh on which all computation are made
  Mesh & mesh;

  /// normals at quadrature points
  ByElementTypeReal normals_on_quad_points;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "fem_inline_impl.cc"
#endif


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const FEM & _this)
{
  _this.printself(stream);
  return stream;
}


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const QuadraturePoint & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#include "fem_template.hh"

#endif /* __AKANTU_FEM_HH__ */
