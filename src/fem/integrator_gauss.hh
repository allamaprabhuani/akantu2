/**
 * @file   integrator_gauss.hh
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Feb 10 16:42:34 2011
 *
 * @brief  Gauss integration facilities
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
#ifndef __AKANTU_INTEGRATOR_GAUSS_HH__
#define __AKANTU_INTEGRATOR_GAUSS_HH__

/* -------------------------------------------------------------------------- */
#include "integrator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class IntegratorGauss : public Integrator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  IntegratorGauss(const Mesh & mesh,
		  const ID & id = "integrator_gauss",
		  const MemoryID & memory_id = 0);

  virtual ~IntegratorGauss(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// precompute jacobians on elements of type "type"
  template <ElementType type>
  void precomputeJacobiansOnQuadraturePoints(const GhostType & ghost_type);


  /// integrate f on the element "elem" of type "type"
  template <ElementType type>
  inline void integrateOnElement(const Vector<Real> & f,
				 Real * intf,
				 UInt nb_degree_of_freedom,
				 const UInt elem,
				 const GhostType & ghost_type) const;

  /// integrate f for all elements of type "type"
  template <ElementType type>
  void integrate(const Vector<Real> & in_f,
		 Vector<Real> &intf,
		 UInt nb_degree_of_freedom,
		 const GhostType & ghost_type,
		 const Vector<UInt> * filter_elements) const;

  /// integrate one element scalar value on all elements of type "type"
  template <ElementType type>
  Real integrate(const types::RVector & in_f,
		 UInt index,
		 const GhostType & ghost_type) const;


  /// integrate scalar field in_f
  template <ElementType type>
  Real integrate(const Vector<Real> & in_f,
		 const GhostType & ghost_type,
		 const Vector<UInt> * filter_elements) const;

  /// integrate partially around a quadrature point (@f$ intf_q = f_q * J_q * w_q @f$)
  template <ElementType type>
  void integrateOnQuadraturePoints(const Vector<Real> & in_f,
				   Vector<Real> &intf,
				   UInt nb_degree_of_freedom,
				   const GhostType & ghost_type,
				   const Vector<UInt> * filter_elements) const;

  /// return a vector with quadrature points natural coordinates
  template <ElementType type>
  const Vector<Real> & getQuadraturePoints(const GhostType & ghost_type) const;

  /// compute the vector of quadrature points natural coordinates
  template <ElementType type> void computeQuadraturePoints(const GhostType & ghost_type);

  /// check that the jacobians are not negative
  template <ElementType type> void checkJacobians(const GhostType & ghost_type) const;

public:

  /// compute the jacobians on quad points for a given element
  template <ElementType type>
  void computeJacobianOnQuadPointsByElement(UInt spatial_dimension,
					    Real * node_coords,
					    UInt nb_nodes_per_element,
					    Real * jacobians);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
private:
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  inline void integrate(Real *f, Real *jac, Real * inte,
			UInt nb_degree_of_freedom,
			UInt nb_quadrature_points) const;

private:

  ByElementTypeReal quadrature_points;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "integrator_gauss_inline_impl.cc"
#endif

__END_AKANTU__

#endif /* __AKANTU_INTEGRATOR_GAUSS_HH__ */
