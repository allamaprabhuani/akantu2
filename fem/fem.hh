/**
 * @file   fem.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jul 16 10:24:24 2010
 *
 * @brief  FEM class
 *
 * @section LICENSE
 *
 * \<insert license here\>
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

  /// pre-compute all the shape functions, their derivatives and the jacobians
  void initShapeFunctions(GhostType ghost_type = _not_ghost);

  /// pre-compute normals on quadrature points
  void computeNormalsOnQuadPoints(GhostType ghost_type = _not_ghost);

  /// interpolate nodal values on the quadrature points
  void interpolateOnQuadraturePoints(const Vector<Real> &u,
				     Vector<Real> &uq,
				     UInt nb_degre_of_freedom,
				     const ElementType & type,
				     GhostType ghost_type = _not_ghost,
				     const Vector<UInt> * filter_elements = NULL) const;

  /// compute the gradient of u on the quadrature points
  void gradientOnQuadraturePoints(const Vector<Real> &u,
				  Vector<Real> &nablauq,
				  UInt nb_degre_of_freedom,
				  const ElementType & type,
				  GhostType ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL) const;

  /// integrate f for all elements of type "type"
  void integrate(const Vector<Real> & f,
		 Vector<Real> &intf,
		 UInt nb_degre_of_freedom,
		 const ElementType & type,
		 GhostType ghost_type = _not_ghost,
		 const Vector<UInt> * filter_elements = NULL) const;

  /// integrate f on the element "elem" of type "type"
  inline void integrate(const Vector<Real> & f,
			Real *intf,
			UInt nb_degre_of_freedom,
			const ElementType & type,
			const UInt elem,
			GhostType ghost_type = _not_ghost) const;

  /// integrate a scalar value on all elements of type "type"
  Real integrate(const Vector<Real> & f,
		 const ElementType & type,
		 GhostType ghost_type = _not_ghost,
		 const Vector<UInt> * filter_elements = NULL) const;


  /// assemble vectors
  void assembleVector(const Vector<Real> & elementary_vect,
		      Vector<Real> & nodal_values,
		      UInt nb_degre_of_freedom,
		      const ElementType & type,
		      GhostType ghost_type = _not_ghost,
		      const Vector<UInt> * filter_elements = NULL,
		      Real scale_factor = 1) const;


  /// assemble matrix in the complete sparse matrix
  void assembleMatrix() {};

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// initialise the class
  void init();

  /// integrate on one element
  inline void integrate(Real *f, Real *jac, Real * inte,
			UInt nb_degre_of_freedom,
			UInt nb_quadrature_points) const;
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(ElementDimension, element_dimension, UInt);

  /// get the mesh contained in the fem object
  inline Mesh & getMesh() const;

  /// get the size of the shapes returned by the element class
  static inline UInt getShapeSize(const ElementType & type);

  /// get the number of quadrature points
  static inline UInt getNbQuadraturePoints(const ElementType & type);

  /// get the size of the shapes derivatives returned by the element class
  static inline UInt getShapeDerivativesSize(const ElementType & type);

  /// get the in-radius of an element
  static inline Real getElementInradius(Real * coord, const ElementType & type);

  /// get a the shapes vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Shapes, shapes, const Vector<Real> &);

  /// get a the shapes derivatives vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ShapesDerivatives, shapes_derivatives, const Vector<Real> &);

  /// get the normals on quadrature points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(NormalsOnQuadPoints, normals_on_quad_points, const Vector<Real> &);

  /// get a the ghost shapes vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostShapes, ghost_shapes,
				   const Vector<Real> &);

  /// get a the ghost shapes derivatives vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostShapesDerivatives, ghost_shapes_derivatives,
				   const Vector<Real> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// id of the fem object
  FEMID id;

  /// spatial dimension of the problem
  UInt element_dimension;

  /// the mesh on which all computation are made
  Mesh * mesh;

  /// shape functions for all elements
  ByElementTypeReal shapes;

  /// shape derivatives for all elements
  ByElementTypeReal shapes_derivatives;

  /// jacobians for all elements
  ByElementTypeReal jacobians;

  /// normals at quadrature points
  ByElementTypeReal normals_on_quad_points;

  /// shape functions for all elements
  ByElementTypeReal ghost_shapes;

  /// shape derivatives for all elements
  ByElementTypeReal ghost_shapes_derivatives;

  /// jacobians for all elements
  ByElementTypeReal ghost_jacobians;
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


#endif /* __AKANTU_FEM_HH__ */
