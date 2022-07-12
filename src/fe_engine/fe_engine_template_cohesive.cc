/**
 * @file   fe_engine_template_cohesive.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Oct 31 2012
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  Specialization of the FEEngineTemplate for cohesive element
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
#include "integrator_gauss.hh"
#include "shape_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include "fe_engine_template.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/* compatibility functions */
/* -------------------------------------------------------------------------- */
template <>
Real FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive,
                      DefaultIntegrationOrderFunctor>::
    integrate(const Array<Real> & f, ElementType type, GhostType ghost_type,
              const Array<Idx> & filter_elements) const {
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
  auto nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  auto nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.size() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID()
                                      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == 1,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
#endif

  Real integral = tuple_dispatch<ElementTypes_t<_ek_cohesive>>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return integrator.integrate<type>(f, ghost_type, filter_elements);
      },
      type);

  AKANTU_DEBUG_OUT();
  return integral;
}

/* -------------------------------------------------------------------------- */
template <>
void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive,
                      DefaultIntegrationOrderFunctor>::
    integrate(const Array<Real> & f, Array<Real> & intf,
              Int nb_degree_of_freedom, ElementType type, GhostType ghost_type,
              const Array<Idx> & filter_elements) const {

#ifndef AKANTU_NDEBUG
  auto nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  auto nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.size() == Int(nb_element * nb_quadrature_points),
                      "The vector f(" << f.getID() << " size " << f.size()
                                      << ") has not the good size ("
                                      << nb_element << ").");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == Int(nb_degree_of_freedom),
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == Int(nb_degree_of_freedom),
                      "The vector intf("
                          << intf.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.size() == nb_element,
                      "The vector intf(" << intf.getID()
                                         << ") has not the good size.");
#endif

  tuple_dispatch<ElementTypes_t<_ek_cohesive>>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        integrator.integrate<type>(f, intf, nb_degree_of_freedom, ghost_type,
                                   filter_elements);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <>
void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive,
                      DefaultIntegrationOrderFunctor>::
    gradientOnIntegrationPoints(const Array<Real> & /* u */,
                                Array<Real> & /*  nablauq */,
                                Int /* nb_degree_of_freedom */,
                                ElementType /* type  */,
                                GhostType /*  ghost_type */,
                                const Array<Idx> &) const {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
