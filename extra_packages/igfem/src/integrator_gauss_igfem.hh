/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_INTEGRATOR_IGFEM_HH_
#define AKANTU_INTEGRATOR_IGFEM_HH_

/* -------------------------------------------------------------------------- */
#include "integrator.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {

template <class IOF> class IntegratorGauss<_ek_igfem, IOF> : public Integrator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  IntegratorGauss(const Mesh & mesh, const ID & id = "integrator_gauss");

  virtual ~IntegratorGauss(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void initIntegrator(const Array<Real> & nodes, ElementType type,
                             GhostType ghost_type);

  /// precompute jacobians on elements of type "type"
  template <ElementType type>
  void precomputeJacobiansOnQuadraturePoints(const Array<Real> & nodes,
                                             GhostType ghost_type);

  /// integrate f on the element "elem" of type "type"
  template <ElementType type>
  inline void integrateOnElement(const Array<Real> & f, Real * intf,
                                 UInt nb_degree_of_freedom, UInt elem,
                                 GhostType ghost_type) const;

  /// integrate f for all elements of type "type"
  template <ElementType type>
  void integrate(const Array<Real> & in_f, Array<Real> & intf,
                 UInt nb_degree_of_freedom, GhostType ghost_type,
                 const Array<Int> & filter_elements) const;

  /// integrate one element scalar value on all elements of type "type"
  template <ElementType type>
  Real integrate(const Vector<Real> & in_f, UInt index,
                 GhostType ghost_type) const;

  /// integrate scalar field in_f
  template <ElementType type>
  Real integrate(const Array<Real> & in_f, GhostType ghost_type,
                 const Array<Int> & filter_elements) const;

  /// integrate partially around a quadrature point (@f$ intf_q = f_q * J_q *
  /// w_q @f$)
  template <ElementType type>
  void
  integrateOnIntegrationPoints(const Array<Real> & in_f, Array<Real> & intf,
                               UInt nb_degree_of_freedom, GhostType ghost_type,
                               const Array<Int> & filter_elements) const;
  /// return a vector with quadrature points natural coordinates
  template <ElementType type>
  const Matrix<Real> & getIntegrationPoints(GhostType ghost_type) const;

  /// return the number of quadrature points for a given element type
  template <ElementType type>
  inline UInt getNbIntegrationPoints(GhostType ghost_type = _not_ghost) const;

  /// compute the vector of quadrature points natural coordinates
  template <ElementType type>
  void computeQuadraturePoints(GhostType ghost_type);

  /// check that the jacobians are not negative
  template <ElementType type> void checkJacobians(GhostType ghost_type) const;

public:
  /// compute the jacobians on quad points for a given element
  template <ElementType type>
  void computeJacobianOnQuadPointsByElement(const Matrix<Real> & node_coords,
                                            Vector<Real> & jacobians);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  inline void integrate(Real * f, Real * jac, Real * inte,
                        UInt nb_degree_of_freedom,
                        UInt nb_quadrature_points) const;

private:
  ElementTypeMap<Matrix<Real>> quadrature_points;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "integrator_gauss_igfem_inline_impl.hh"

} // namespace akantu

#endif /*AKANTU_INTEGRATOR_IGFEM_HH_ */
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
