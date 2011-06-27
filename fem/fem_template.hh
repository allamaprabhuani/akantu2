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
	      FEMID id = "fem", MemoryID memory_id = 0);

  virtual ~FEMTemplate();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// pre-compute all the shape functions, their derivatives and the jacobians
  void initShapeFunctions(GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */

  /// integrate f for all elements of type "type"
  void integrate(const Vector<Real> & f,
		 Vector<Real> &intf,
		 UInt nb_degre_of_freedom,
   		 const ElementType & type,
   		 GhostType ghost_type = _not_ghost,
   		 const Vector<UInt> * filter_elements = NULL) const;

  /// integrate a scalar value on all elements of type "type"
  Real integrate(const Vector<Real> & f,
   		 const ElementType & type,
   		 GhostType ghost_type = _not_ghost,
   		 const Vector<UInt> * filter_elements = NULL) const;

  /// get the number of quadrature points
  UInt getNbQuadraturePoints(const ElementType & type);

  /// get shapes precomputed
  const Vector<Real> & getShapes(const ElementType & type, 
				 const GhostType & ghost_type);
  /// get the derivatives of shapes
  const Vector<Real> & getShapesDerivatives(const ElementType & type,
					    const GhostType & ghost_type);
  /// get quadrature points
  const Vector<Real> & getQuadraturePoints(const ElementType & type);


  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */


  void gradientOnQuadraturePoints(const Vector<Real> &u,
				  Vector<Real> &nablauq,
				  const UInt nb_degre_of_freedom,
				  const ElementType & type,
				  GhostType ghost_type = _not_ghost,
				  const Vector<UInt> * filter_elements = NULL);

  void interpolateOnQuadraturePoints(const Vector<Real> &u,
				     Vector<Real> &uq,
				     UInt nb_degre_of_freedom,
				     const ElementType & type,
				     GhostType ghost_type = _not_ghost,
				     const Vector<UInt> * filter_elements = NULL) const;

  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */

  /// pre-compute normals on control points
  void computeNormalsOnControlPoints(GhostType ghost_type = _not_ghost);

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const{};

  void assembleFieldLumped(const Vector<Real> & field_1,
			   UInt nb_degree_of_freedom,
			   Vector<Real> & lumped,
			   const Vector<Int> & equation_number,
			   ElementType type,
			   GhostType ghost_type);

  void assembleFieldMatrix(const Vector<Real> & field,
			   UInt nb_degree_of_freedom,
			   SparseMatrix & matrix,
			   ElementType type,
			   GhostType ghost_type);


private:

  template <ElementType type>
  void assembleLumpedTemplate(const Vector<Real> & field_1,
			      UInt nb_degree_of_freedom,
			      Vector<Real> & lumped,
			      const Vector<Int> & equation_number,
			      GhostType ghost_type);
  template <ElementType type>
  void assembleLumpedRowSum(const Vector<Real> & field_1,
			    UInt nb_degree_of_freedom,
			    Vector<Real> & lumped,
			    const Vector<Int> & equation_number,
			    GhostType ghost_type);
  template <ElementType type>
  void assembleLumpedDiagonalScaling(const Vector<Real> & field_1,
				     UInt nb_degree_of_freedom,
				     Vector<Real> & lumped,
				     const Vector<Int> & equation_number,				
				     GhostType ghost_type);

  template <ElementType type>
  void assembleFieldMatrix(const Vector<Real> & field,
			   UInt nb_degree_of_freedom,
			   SparseMatrix & matrix,
			   GhostType ghost_type);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

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

//#include "fem_template_inline_impl.cc"

/// standard output stream operator

// inline std::ostream & operator <<(std::ostream & stream, const FEMTemplate & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__


#endif /* __AKANTU_FEM_TEMPLATE_HH__ */
