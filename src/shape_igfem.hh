/**
 * @file   shape_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  shape functions for interface-enriched generalized FEM
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
//#include "aka_array.hh"
#include "shape_functions.hh"

#ifndef __AKANTU_SHAPE_IGFEM_HH__
#define __AKANTU_SHAPE_IGFEM_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <>
class ShapeLagrange<_ek_igfem> : public ShapeFunctions{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ShapeLagrange(const Mesh & mesh,
		const ID & id = "shape_igfem",
		const MemoryID & memory_id = 0);
  

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
inline void initShapeFunctions(const Array<Real> & nodes,
				 const Matrix<Real> & control_points,
				 const Matrix<Real> & control_points_1,
				 const Matrix<Real> & control_points_2,
				 const ElementType & type,
				 const GhostType & ghost_type);


  /// pre compute all shapes on the element control points from natural coordinates
  template<ElementType type>
  void precomputeShapesOnControlPoints(const Array<Real> & nodes,
				       GhostType ghost_type);

  /// pre compute all shape derivatives on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnControlPoints(const Array<Real> & nodes,
						 GhostType ghost_type);

  /// interpolate nodal values on the control points
  template <ElementType type>
  void interpolateOnControlPoints(const Array<Real> &u,
				  Array<Real> &uq,
				  UInt nb_degree_of_freedom,
				  GhostType ghost_type = _not_ghost,
				  const Array<UInt> & filter_elements = empty_filter) const;

  /// compute the gradient of u on the control points
  template <ElementType type>
  void gradientOnControlPoints(const Array<Real> &u,
			       Array<Real> &nablauq,
			       UInt nb_degree_of_freedom,
			       GhostType ghost_type = _not_ghost,
			       const Array<UInt> & filter_elements = empty_filter) const;

  /// multiply a field by shape functions  @f$ fts_{ij} = f_i * \varphi_j @f$
  template <ElementType type>
  void fieldTimesShapes(const Array<Real> & field,
			Array<Real> & field_times_shapes,
			GhostType ghost_type) const;

  /// find natural coords from real coords provided an element
  template <ElementType type>
  void inverseMap(const Vector<Real> & real_coords,
		  UInt element,
		  Vector<Real> & natural_coords,
		  const GhostType & ghost_type = _not_ghost) const;

  /// return true if the coordinates provided are inside the element, false otherwise
  template <ElementType type>
  bool contains(const Vector<Real> & real_coords,
		UInt elem,
		const GhostType & ghost_type) const;

  /// compute the shape on a provided point
  template <ElementType type>
  void computeShapes(const Vector<Real> & real_coords,
		     UInt elem,
		     Vector<Real> & shapes,
		     const GhostType & ghost_type) const;

  /// compute the shape derivatives on a provided point
  template <ElementType type>
  void computeShapeDerivatives(const Matrix<Real> & real_coords,
		     UInt elem,
		     Tensor3<Real> & shapes,
		     const GhostType & ghost_type) const;

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// compute the shape derivatives on control points for a given element
  template <ElementType type>
  inline void computeShapeDerivativesOnCPointsByElement(const Matrix<Real> & node_coords,
							const Matrix<Real> & natural_coords,
							Tensor3<Real> & shapesd) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get a the shapes vector
  inline const Array<Real> & getShapes(const ElementType & el_type,
					const GhostType & ghost_type = _not_ghost) const;

  /// get a the shapes derivatives vector
  inline const Array<Real> & getShapesDerivatives(const ElementType & el_type,
						   const GhostType & ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// shape functions for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes;

  /// shape functions derivatives for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes_derivatives;

  /// additional control points for the IGFEM formulation
  ElementTypeMapArray<Real> igfem_control_points;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "shape_igfem_inline_impl.cc"
/// standard output stream operator
// template <class ShapeFunction>
// inline std::ostream & operator <<(std::ostream & stream, const ShapeIGFEM<ShapeFunction> & _this)
// {
//   _this.printself(stream);
//   return stream;
// }
__END_AKANTU__

#endif /* __AKANTU_SHAPE_IGFEM_HH__ */
