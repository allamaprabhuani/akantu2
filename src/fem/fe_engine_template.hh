/**
 * @file   fe_engine_template.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date   Tue Feb 15 16:32:44 2011
 *
 * @brief  templated class that calls integration and shape objects
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

#ifndef __AKANTU_FE_ENGINE_TEMPLATE_HH__
#define __AKANTU_FE_ENGINE_TEMPLATE_HH__

/* -------------------------------------------------------------------------- */
#include "fe_engine.hh"
#include "integrator.hh"
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__
template<ElementKind k> class AssembleLumpedTemplateHelper;


/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind = _ek_regular>
class FEEngineTemplate : public FEEngine {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef I<kind> Integ;
  typedef S<kind> Shape;

  FEEngineTemplate(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
		   ID id = "fem", MemoryID memory_id = 0);

  virtual ~FEEngineTemplate();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// pre-compute all the shape functions, their derivatives and the jacobians
  void initShapeFunctions(const GhostType & ghost_type = _not_ghost);
  void initShapeFunctions(const Array<Real> & nodes,
			  const GhostType & ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */
  /// integrate f for all elements of type "type"
  void integrate(const Array<Real> & f,
		 Array<Real> &intf,
		 UInt nb_degree_of_freedom,
   		 const ElementType & type,
   		 const GhostType & ghost_type = _not_ghost,
   		 const Array<UInt> & filter_elements = empty_filter) const;

  /// integrate a scalar value on all elements of type "type"
  Real integrate(const Array<Real> & f,
   		 const ElementType & type,
   		 const GhostType & ghost_type = _not_ghost,
   		 const Array<UInt> & filter_elements = empty_filter) const;

  /// integrate one element scalar value on all elements of type "type"
  virtual Real integrate(const Vector<Real> & f,
			 const ElementType & type,
			 UInt index, const GhostType & ghost_type = _not_ghost) const;


  /// integrate partially around a quadrature point (@f$ intf_q = f_q * J_q * w_q @f$)
  void integrateOnQuadraturePoints(const Array<Real> & f,
				   Array<Real> &intf,
				   UInt nb_degree_of_freedom,
				   const ElementType & type,
				   const GhostType & ghost_type = _not_ghost,
				   const Array<UInt> & filter_elements = empty_filter) const;


  /// get the number of quadrature points
  UInt getNbQuadraturePoints(const ElementType & type,
			     const GhostType & ghost_type = _not_ghost) const;

  /// get shapes precomputed
  const Array<Real> & getShapes(const ElementType & type,
				const GhostType & ghost_type = _not_ghost) const;

  /// get the derivatives of shapes
  const Array<Real> & getShapesDerivatives(const ElementType & type,
					   const GhostType & ghost_type = _not_ghost,
					   UInt id=0) const;

  /// get quadrature points
  const inline Matrix<Real> & getQuadraturePoints(const ElementType & type,
						  const GhostType & ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */

  /// compute the gradient of a nodal field on the quadrature points
  void gradientOnQuadraturePoints(const Array<Real> &u,
				  Array<Real> &nablauq,
				  const UInt nb_degree_of_freedom,
				  const ElementType & type,
				  const GhostType & ghost_type = _not_ghost,
				  const Array<UInt> & filter_elements = empty_filter) const;

  /// interpolate a nodal field on the quadrature points
  void interpolateOnQuadraturePoints(const Array<Real> &u,
				     Array<Real> &uq,
				     UInt nb_degree_of_freedom,
				     const ElementType & type,
				     const GhostType & ghost_type = _not_ghost,
				     const Array<UInt> & filter_elements = empty_filter) const;

  /// interpolate a nodal field on the quadrature points given a by_element_type
  void interpolateOnQuadraturePoints(const Array<Real> & u,
				     ElementTypeMapArray<Real> & uq,
                                     const ElementTypeMapArray<UInt> * filter_elements = NULL) const;

  /// find natural coords from real coords provided an element
  void inverseMap(const Vector<Real> & real_coords,
		  UInt element,
		  const ElementType & type,
		  Vector<Real> & natural_coords,
		  const GhostType & ghost_type = _not_ghost) const;

  /// return true if the coordinates provided are inside the element, false otherwise
  inline bool contains(const Vector<Real> & real_coords,
		       UInt element,
		       const ElementType & type,
		       const GhostType & ghost_type = _not_ghost) const;

  /// compute the shape on a provided point
  inline void computeShapes(const Vector<Real> & real_coords,
			    UInt element,
			    const ElementType & type,
			    Vector<Real> & shapes,
			    const GhostType & ghost_type = _not_ghost) const;



  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */

  /// pre-compute normals on control points
  void computeNormalsOnControlPoints(const GhostType & ghost_type = _not_ghost);
  void computeNormalsOnControlPoints(const Array<Real> & field,
				     const GhostType & ghost_type = _not_ghost);
  void computeNormalsOnControlPoints(const Array<Real> & field,
				     Array<Real> & normal,
				     const ElementType & type,
				     const GhostType & ghost_type = _not_ghost) const;
  template<ElementType type>
  void computeNormalsOnControlPoints(const Array<Real> & field,
				     Array<Real> & normal,
				     const GhostType & ghost_type) const;


  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const{};

  void assembleFieldLumped(const Array<Real> & field_1,
			   UInt nb_degree_of_freedom,
			   Array<Real> & lumped,
			   const Array<Int> & equation_number,
			   ElementType type,
			   const GhostType & ghost_type = _not_ghost) const;

  void assembleFieldMatrix(const Array<Real> & field,
			   UInt nb_degree_of_freedom,
			   SparseMatrix & matrix,
			   ElementType type,
			   const GhostType & ghost_type = _not_ghost) const;

private:
  friend class AssembleLumpedTemplateHelper<kind>;


  template <ElementType type>
  void assembleLumpedTemplate(const Array<Real> & field_1,
			      UInt nb_degree_of_freedom,
			      Array<Real> & lumped,
			      const Array<Int> & equation_number,
			      const GhostType & ghost_type) const;

  /// @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV = \int \rho \varphi_i dV @f$
  template <ElementType type>
  void assembleLumpedRowSum(const Array<Real> & field_1,
			    UInt nb_degree_of_freedom,
			    Array<Real> & lumped,
			    const Array<Int> & equation_number,
			    const GhostType & ghost_type) const;

  /// @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
  template <ElementType type>
  void assembleLumpedDiagonalScaling(const Array<Real> & field_1,
				     UInt nb_degree_of_freedom,
				     Array<Real> & lumped,
				     const Array<Int> & equation_number,
				     const GhostType & ghost_type) const;

  template <ElementType type>
  void assembleFieldMatrix(const Array<Real> & field,
			   UInt nb_degree_of_freedom,
			   SparseMatrix & matrix,
			   const GhostType & ghost_type) const;


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  const ShapeFunctions & getShapeFunctionsInterface() const { return shape_functions; };
  const Shape & getShapeFunctions() const { return shape_functions; };

  const Integrator & getIntegratorInterface() const { return integrator; };
  const Integ & getIntegrator() const { return integrator; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Integ integrator;
  Shape shape_functions;
};


__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "fe_engine_template_tmpl.hh"

#endif /* __AKANTU_FE_ENGINE_TEMPLATE_HH__ */
