/**
 * @file   integrator_gauss.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  Gauss integration facilities
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "integrator.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_INTEGRATOR_GAUSS_HH_
#define AKANTU_INTEGRATOR_GAUSS_HH_

namespace akantu {
namespace integrator {
  namespace details {
    template <ElementKind> struct GaussIntegratorComputeJacobiansHelper;
  } // namespace details
} // namespace integrator

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss : public Integrator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  IntegratorGauss(const Mesh & mesh, Int spatial_dimension,
                  const ID & id = "integrator_gauss");

  ~IntegratorGauss() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initIntegrator(const Array<Real> & nodes, ElementType type,
                      GhostType ghost_type);

  template <ElementType type>
  inline void initIntegrator(const Array<Real> & nodes, GhostType ghost_type);

  /// integrate f on the element "elem" of type "type"
  template <ElementType type>
  inline void integrateOnElement(const Array<Real> & f, Real * intf,
                                 Int nb_degree_of_freedom, const Idx elem,
                                 GhostType ghost_type) const;

  /// integrate f for all elements of type "type"
  template <ElementType type>
  void integrate(const Array<Real> & in_f, Array<Real> & intf,
                 Int nb_degree_of_freedom, GhostType ghost_type,
                 const Array<Idx> & filter_elements) const;

  /// integrate scalar field in_f
  template <ElementType type, Int polynomial_degree>
  Real integrate(const Array<Real> & in_f,
                 GhostType ghost_type = _not_ghost) const;

  /// integrate partially around a quadrature point (@f$ intf_q = f_q * J_q *
  /// w_q @f$)
  template <ElementType type>
  Real integrate(const Vector<Real> & in_f, Idx index,
                 GhostType ghost_type) const;

  /// integrate scalar field in_f
  template <ElementType type>
  Real integrate(const Array<Real> & in_f, GhostType ghost_type,
                 const Array<Idx> & filter_elements) const;

  /// integrate a field without using the pre-computed values
  template <ElementType type, Int polynomial_degree>
  void integrate(const Array<Real> & in_f, Array<Real> & intf,
                 Int nb_degree_of_freedom, GhostType ghost_type) const;

  /// integrate partially around a quadrature point (@f$ intf_q = f_q * J_q *
  /// w_q @f$)
  template <ElementType type>
  void
  integrateOnIntegrationPoints(const Array<Real> & in_f, Array<Real> & intf,
                               Int nb_degree_of_freedom, GhostType ghost_type,
                               const Array<Idx> & filter_elements) const;

  /// return a matrix with quadrature points natural coordinates
  template <ElementType type>
  const Matrix<Real> & getIntegrationPoints(GhostType ghost_type) const;

  /// return number of quadrature points
  template <ElementType type>
  Int getNbIntegrationPoints(GhostType ghost_type) const;

  template <ElementType type, Int n> Matrix<Real> getIntegrationPoints() const;
  template <ElementType type, Int n> Vector<Real> getIntegrationWeights() const;

protected:
  friend struct integrator::details::GaussIntegratorComputeJacobiansHelper<
      kind>;

  template <ElementType type>
  void computeJacobiansOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & quad_points,
      Array<Real> & jacobians, GhostType ghost_type,
      const Array<Idx> & filter_elements = empty_filter) const;

  void computeJacobiansOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & quad_points,
      Array<Real> & jacobians, ElementType type, GhostType ghost_type,
      const Array<Idx> & filter_elements = empty_filter) const;

  /// precompute jacobians on elements of type "type"
  template <ElementType type>
  void precomputeJacobiansOnQuadraturePoints(const Array<Real> & nodes,
                                             GhostType ghost_type);

  // multiply the jacobians by the integration weights and stores the results in
  // jacobians
  template <ElementType type, Int polynomial_degree>
  void multiplyJacobiansByWeights(
      Array<Real> & jacobians,
      const Array<Idx> & filter_elements = empty_filter) const;

  /// compute the vector of quadrature points natural coordinates
  template <ElementType type>
  void computeQuadraturePoints(GhostType ghost_type);

  /// check that the jacobians are not negative
  template <ElementType type> void checkJacobians(GhostType ghost_type) const;

  /// internal integrate partially around a quadrature point (@f$ intf_q = f_q *
  /// J_q *
  /// w_q @f$)
  void integrateOnIntegrationPoints(const Array<Real> & in_f,
                                    Array<Real> & intf,
                                    Int nb_degree_of_freedom,
                                    const Array<Real> & jacobians,
                                    Int nb_element) const;

  void integrate(const Array<Real> & in_f, Array<Real> & intf,
                 Int nb_degree_of_freedom, const Array<Real> & jacobians,
                 Int nb_element) const;

public:
  /// compute the jacobians on quad points for a given element
  template <ElementType type, class D1, class D2, class D3>
  inline void computeJacobianOnQuadPointsByElement(
      const Eigen::MatrixBase<D1> & node_coords,
      const Eigen::MatrixBase<D2> & quad,
      Eigen::MatrixBase<D3> & jacobians) const;

public:
  void onElementsAdded(const Array<Element> & elements) override;

  template <ElementType type>
  void onElementsAddedByType(const Array<Idx> & new_elements,
                             GhostType ghost_type);

  decltype(auto) getIntegrationPoints() { return (quadrature_points); }
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// integrate the field f with the jacobian jac -> inte
  inline void integrate(const Real * f, const Real * jac, Real * inte,
                        Int nb_degree_of_freedom,
                        Int nb_quadrature_points) const;

private:
  /// ElementTypeMap of the quadrature points
  ElementTypeMap<Matrix<Real>> quadrature_points;
};

} // namespace akantu

#include "integrator_gauss_inline_impl.hh"

#endif /* AKANTU_INTEGRATOR_GAUSS_HH_ */
