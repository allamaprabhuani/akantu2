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

#ifndef __MYFEM_FEM_HH__
#define __MYFEM_FEM_HH__

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "memory.hh"
#include "mesh.hh"
#include "element_class.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef Vector<Real> * ConnectivityTypeDataReal[_max_element_type];

class FEM : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FEM(UInt spatial_dimension, FEMID id = "fem", MemoryID memory_id = 0);

  FEM(Mesh & mesh, UInt spatial_dimension = 0,
      FEMID id = "fem", MemoryID memory_id = 0);

  virtual ~FEM();

  typedef std::map<ElementType, Vector<Real> *> ByTypeRealMap;
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// pre-compute all the shape functions, their derivatives and the jacobians
  void initShapeFunctions();

  /// compute the  volume of an element
  Real volume(ElementType type, Int element);

  /// interpolate nodal values on quadrature points
  void interpolateOnQuadraturePoints(const Vector<Real> &inval,
				     Vector<Real> &valonquad,
				     ElementType type,
				     const Vector<UInt> * element = NULL);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  inline Mesh & getMesh() const;

  /// get the number of quadrature points of an element
  inline UInt getNbQuadraturePoints(ElementType type) const;

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
  ConnectivityTypeDataReal shapes;

  /// shape derivatives for all elements
  ConnectivityTypeDataReal shapes_derivatives;

  /// jacobians for all elements
  ConnectivityTypeDataReal jacobians;

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


#endif /* __MYFEM_FEM_HH__ */
