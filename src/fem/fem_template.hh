/**
 * @file   fem_template.hh
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Thu Feb 10 10:55:21 2011
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

#ifndef __AKANTU_FEM_TEMPLATE_HH__
#define __AKANTU_FEM_TEMPLATE_HH__

/* -------------------------------------------------------------------------- */
#include "fem.hh"
#include "integrator.hh"
#include "shape_functions.hh"
#include "shape_lagrange.hh"
#include "shape_cohesive.hh"
#include "shape_linked.hh"
#include "integrator_gauss.hh"
#include "integrator_cohesive.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <typename Integ, typename Shape>
class FEMTemplate : public FEM{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FEMTemplate(Mesh & mesh, UInt spatial_dimension = 0,
	      ID id = "fem", MemoryID memory_id = 0);

  virtual ~FEMTemplate();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// pre-compute all the shape functions, their derivatives and the jacobians
  void initShapeFunctions(const GhostType & ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */

  /// integrate f for all elements of type "type"
  void integrate(const Vector<Real> & f,
		 Vector<Real> &intf,
		 UInt nb_degree_of_freedom,
   		 const ElementType & type,
   		 const GhostType & ghost_type = _not_ghost,
   		 const Vector<UInt> * filter_elements = NULL) const;

  /// integrate a scalar value on all elements of type "type"
  Real integrate(const Vector<Real> & f,
   		 const ElementType & type,
   		 const GhostType & ghost_type = _not_ghost,
   		 const Vector<UInt> * filter_elements = NULL) const;

  /// integrate partially around a quadrature point (@f$ intf_q = f_q * J_q * w_q @f$)
  void integrateOnQuadraturePoints(const Vector<Real> & f,
				   Vector<Real> &intf,
				   UInt nb_degree_of_freedom,
				   const ElementType & type,
				   const GhostType & ghost_type = _not_ghost,
				   const Vector<UInt> * filter_elements = NULL) const;


  /// get the number of quadrature points
  UInt getNbQuadraturePoints(const ElementType & type,
			     const GhostType & ghost_type = _not_ghost) const;

  /// get shapes precomputed
  const Vector<Real> & getShapes(const ElementType & type,
				 const GhostType & ghost_type = _not_ghost) const;

  /// get the derivatives of shapes
  const Vector<Real> & getShapesDerivatives(const ElementType & type,
					    const GhostType & ghost_type = _not_ghost,
					    UInt id=0) const;

  /// get quadrature points
  const Vector<Real> & getQuadraturePoints(const ElementType & type,
					   const GhostType & ghost_type = _not_ghost) const;


  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */


  void gradientOnQuadraturePoints(const Vector<Real> &u,
				  Vector<Real> &nablauq,
				  const UInt nb_degree_of_freedom,
				  const ElementType & type,
				  const GhostType & ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL) const;

  void interpolateOnQuadraturePoints(const Vector<Real> &u,
				     Vector<Real> &uq,
				     UInt nb_degree_of_freedom,
				     const ElementType & type,
				     const GhostType & ghost_type = _not_ghost,
				     const Vector<UInt> * filter_elements = NULL) const;

  /// find natural coords from real coords provided an element
  void inverseMap(const types::RVector & real_coords,
		  UInt element,
		  const ElementType & type,
		  types::RVector & natural_coords,
		  const GhostType & ghost_type = _not_ghost) const;

  /// return true if the coordinates provided are inside the element, false otherwise
  inline bool contains(const types::RVector & real_coords,
		       UInt element,
		       const ElementType & type,
		       const GhostType & ghost_type = _not_ghost) const;

  /// compute the shape on a provided point
  inline void computeShapes(const types::RVector & real_coords,
			    UInt element,
			    const ElementType & type,
			    types::RVector & shapes,
			    const GhostType & ghost_type = _not_ghost) const;



  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */

  /// pre-compute normals on control points
  void computeNormalsOnControlPoints(const GhostType & ghost_type = _not_ghost);
  void computeNormalsOnControlPoints(const Vector<Real> & field,
				     const GhostType & ghost_type = _not_ghost);
  void computeNormalsOnControlPoints(const Vector<Real> & field,
				     Vector<Real> & normal,
				     const ElementType & type,
				     const GhostType & ghost_type = _not_ghost) const;


  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const{};

  void assembleFieldLumped(const Vector<Real> & field_1,
			   UInt nb_degree_of_freedom,
			   Vector<Real> & lumped,
			   const Vector<Int> & equation_number,
			   ElementType type,
			   const GhostType & ghost_type = _not_ghost) const;

  void assembleFieldMatrix(const Vector<Real> & field,
			   UInt nb_degree_of_freedom,
			   SparseMatrix & matrix,
			   ElementType type,
			   const GhostType & ghost_type = _not_ghost) const;


private:

  template <ElementType type>
  void assembleLumpedTemplate(const Vector<Real> & field_1,
			      UInt nb_degree_of_freedom,
			      Vector<Real> & lumped,
			      const Vector<Int> & equation_number,
			      const GhostType & ghost_type) const;

  /// @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV = \int \rho \varphi_i dV @f$
  template <ElementType type>
  void assembleLumpedRowSum(const Vector<Real> & field_1,
			    UInt nb_degree_of_freedom,
			    Vector<Real> & lumped,
			    const Vector<Int> & equation_number,
			    const GhostType & ghost_type) const;

  /// @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
  template <ElementType type>
  void assembleLumpedDiagonalScaling(const Vector<Real> & field_1,
				     UInt nb_degree_of_freedom,
				     Vector<Real> & lumped,
				     const Vector<Int> & equation_number,
				     const GhostType & ghost_type) const;

  template <ElementType type>
  void assembleFieldMatrix(const Vector<Real> & field,
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


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "fem_template_inline_impl.cc"

/// standard output stream operator

// inline std::ostream & operator <<(std::ostream & stream, const FEMTemplate & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__


#endif /* __AKANTU_FEM_TEMPLATE_HH__ */
