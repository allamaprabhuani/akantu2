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
#include "common.hh"
#include "memory.hh"
#include "mesh.hh"
#include "element_class.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef Vector<Real> * ByConnectivityTypeReal[_max_element_type];

typedef Vector<Int>  * ByConnectivityTypeInt[_max_element_type];


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
  void initShapeFunctions();

  /// compute the  volume of an element
  Real volume(ElementType type, Int element);

  /// interpolate nodal values on the quadrature points
  void interpolateOnQuadraturePoints(const Vector<Real> &u,
				     Vector<Real> &uq,
				     UInt nb_degre_of_freedom,
				     const ElementType & type,
				     const Vector<UInt> * element = NULL);


  /// compute the gradient of u on the quadrature points
  void gradientOnQuadraturePoints(const Vector<Real> &u,
				  Vector<Real> &nablauq,
				  UInt nb_degre_of_freedom,
				  const ElementType & type,
				  const Vector<UInt> * element = NULL);

  /// integrate f on all elements of type "type"
  void integrate(const Vector<Real> & f,
		 Vector<Real> &intf,
		 UInt nb_degre_of_freedom,
		 const ElementType & type,
		 const Vector<UInt> * filter_elements = NULL);

  /// integrate f on the element "elem" of type "type"
  inline void integrate(const Vector<Real> & f,
			Real *intf,
			UInt nb_degre_of_freedom,
			const ElementType & type,
			const UInt elem);

  /// assemble vectors
  void assembleVector(const Vector<Real> & elementary_vect,
		      Vector<Real> & nodal_values,
		      UInt nb_degre_of_freedom,
		      const ElementType & type,
		      const Vector<UInt> * filter_elements = NULL);

  void assembleMatrix() {};

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  inline void integrate(Real *f, Real *jac, Real * inte,
			UInt nb_degre_of_freedom,
			UInt nb_quadrature_points);
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the mesh contained in the fem object
  inline Mesh & getMesh() const;

  /// get the types list containted in the mesh
  inline const Mesh::ConnectivityTypeList & getConnectivityTypeList() const;

  /// get the number of nodes in the mesh
  inline UInt getNbNodes() const;

  /// get the number of element of a type in the mesh
  inline UInt getNbElement(const ElementType & type) const;

  /// get the number of quadrature points of an element
  inline UInt getNbQuadraturePoints(const ElementType & type) const;

  /// get spatial dimension of a type of element
  inline UInt getSpatialDimension(const ElementType & type) const;

  /// get the number of nodes per element for a given element type
  inline UInt getNbNodesPerElement(const ElementType & type) const;

  /// get a the shape vector
  inline const Vector<Real> & getShapes(const ElementType & type) const;

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
  ByConnectivityTypeReal shapes;

  /// shape derivatives for all elements
  ByConnectivityTypeReal shapes_derivatives;

  /// jacobians for all elements
  ByConnectivityTypeReal jacobians;

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
