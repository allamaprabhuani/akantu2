/**
 * @file   fem.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jul 16 10:24:24 2010
 *
 * @brief  FEM class
 *
 * @section LICENSE
 *
 * <insert license here>
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

typedef Vector<Real> * ByElementTypeReal[_max_element_type];
typedef Vector<Int>  * ByElementTypeInt[_max_element_type];
typedef Vector<UInt> * ByElementTypeUInt[_max_element_type];

class FEM : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FEM(UInt spatial_dimension, FEMID id = "fem", MemoryID memory_id = 0);

  FEM(Mesh & mesh, UInt spatial_dimension = 0,
      FEMID id = "fem", MemoryID memory_id = 0);

  virtual ~FEM();

  //  typedef std::map<ElementType, Vector<Real> *> ByTypeRealMap;
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// pre-compute all the shape functions, their derivatives and the jacobians
  void initShapeFunctions(bool local = true);

  /// compute the  volume of an element
  Real volume(ElementType type, Int element) const;

  /// interpolate nodal values on the quadrature points
  void interpolateOnQuadraturePoints(const Vector<Real> &u,
				     Vector<Real> &uq,
				     UInt nb_degre_of_freedom,
				     const ElementType & type,
				     bool local = true,
				     const Vector<UInt> * filter_elements = NULL) const;


  /// compute the gradient of u on the quadrature points
  void gradientOnQuadraturePoints(const Vector<Real> &u,
				  Vector<Real> &nablauq,
				  UInt nb_degre_of_freedom,
				  const ElementType & type,
				  bool local = true,
				  const Vector<UInt> * filter_elements = NULL) const;

  /// integrate f for all elements of type "type"
  void integrate(const Vector<Real> & f,
		 Vector<Real> &intf,
		 UInt nb_degre_of_freedom,
		 const ElementType & type,
		 bool local = true,
		 const Vector<UInt> * filter_elements = NULL) const;

  /// integrate f on the element "elem" of type "type"
  inline void integrate(const Vector<Real> & f,
			Real *intf,
			UInt nb_degre_of_freedom,
			const ElementType & type,
			const UInt elem,
			bool local = true) const;

  /// integrate a scalar value on all elements of type "type"
  Real integrate(const Vector<Real> & f,
		 const ElementType & type,
		 bool local = true,
		 const Vector<UInt> * filter_elements = NULL) const;


  /// assemble vectors
  void assembleVector(const Vector<Real> & elementary_vect,
		      Vector<Real> & nodal_values,
		      UInt nb_degre_of_freedom,
		      const ElementType & type,
		      bool local = true,
		      const Vector<UInt> * filter_elements = NULL,
		      Real scale_factor = 1) const;

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

  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the mesh contained in the fem object
  inline Mesh & getMesh() const;

  /// get the number of quadrature points of an element
  static inline UInt getNbQuadraturePoints(const ElementType & type);

  /// get the size of the shapes returned by the element class
  static inline UInt getShapeSize(const ElementType & type);

  /// get the size of the shapes derivatives returned by the element class
  static inline UInt getShapeDerivativesSize(const ElementType & type);

  /// get the size of the jacobian returned by the element class
  static inline UInt getJacobianSize(const ElementType & type);

  /// get the in-radius of an element
  static inline Real getElementInradius(Real * coord, const ElementType & type);

  /// get a the shapes vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Shapes, shapes, const Vector<Real> &);

  /// get a the shapes derivatives vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ShapesDerivatives, shapes_derivatives, const Vector<Real> &);

#ifdef AKANTU_USE_MPI
  /// get a the ghost shapes vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostShapes, ghost_shapes,
				   const Vector<Real> &);

  /// get a the ghost shapes derivatives vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostShapesDerivatives, ghost_shapes_derivatives,
				   const Vector<Real> &);
#endif //AKANTU_USE_MPI

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// id of the fem object
  FEMID id;

  /// spatial dimension of the problem
  UInt spatial_dimension;

  /// the mesh on which all computation are made
  Mesh * mesh;

  /// has the mesh been created by this object
  bool created_mesh;

  /// shape functions for all elements
  ByElementTypeReal shapes;

  /// shape derivatives for all elements
  ByElementTypeReal shapes_derivatives;

  /// jacobians for all elements
  ByElementTypeReal jacobians;

#ifdef AKANTU_USE_MPI
  /// shape functions for all elements
  ByElementTypeReal ghost_shapes;

  /// shape derivatives for all elements
  ByElementTypeReal ghost_shapes_derivatives;

  /// jacobians for all elements
  ByElementTypeReal ghost_jacobians;
#endif //AKANTU_USE_MPI

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
