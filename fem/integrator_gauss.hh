/**
 * @file   integrator_gauss.hh
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
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

/* -------------------------------------------------------------------------- */
#include "integrator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class IntegratorGauss : public Integrator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  IntegratorGauss(Mesh & mesh, IntegratorID id="IntegratorGauss");
  virtual ~IntegratorGauss(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// precompute jacobians on elements of type "type"
  template <ElementType type>
  void precomputeJacobiansOnQuadraturePoints(const UInt dimension,
					     GhostType ghost_type);


  /// integrate f on the element "elem" of type "type"
  template <ElementType type>
  inline void integrateOnElement(const Vector<Real> & f,
				 Real * intf,
				 UInt nb_degre_of_freedom,
				 const UInt elem,
				 GhostType ghost_type) const;

  /// integrate f for all elements of type "type"
  template <ElementType type>
  void integrate(const Vector<Real> & in_f,
		 Vector<Real> &intf,
		 UInt nb_degre_of_freedom,
		 GhostType ghost_type,
		 const Vector<UInt> * filter_elements) const;

  /// integrate scalar field in_f
  template <ElementType type>
  Real integrate(const Vector<Real> & in_f,
		 GhostType ghost_type,
		 const Vector<UInt> * filter_elements) const;

  /// return a vector with quadrature points natural coordinates
  template <ElementType type> Vector<Real> & getQuadraturePoints() const;

  /// compute the vector of quadrature points natural coordinates
  template <ElementType type> void computeQuadraturePoints();

  // /// function to print the contain of the class
  // virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  // /// get the number of quadrature points
  // static inline UInt getNbQuadraturePoints(const ElementType & type);

public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  inline void integrate(Real *f, Real *jac, Real * inte,
			UInt nb_degre_of_freedom,
			UInt nb_quadrature_points) const;


  ByElementTypeReal quadrature_points;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "integrator_gauss_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const IntegratorGauss & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__
