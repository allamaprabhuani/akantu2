/**
 * @file   fem_template.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Feb 15 16:32:44 2011
 *
 * @brief  implementation of the generic FEMTemplate class
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
#include "integrator_gauss.hh"
#include "fem.hh"
#include "aka_common.hh"
#include "shape_lagrange.hh"
#include "shape_linked.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
//template <typename Integ, typename Shape>
template <>
template <>
void FEMTemplate<IntegratorGauss,ShapeLagrange>::
assembleLumpedTemplate<_triangle_6>(const Vector<Real> & field_1,
				    UInt nb_degree_of_freedom,
				    Vector<Real> & lumped,
				    const Vector<Int> & equation_number,
				    const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_triangle_6>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
void FEMTemplate<IntegratorGauss,ShapeLagrange>::
assembleLumpedTemplate<_tetrahedron_10>(const Vector<Real> & field_1,
					UInt nb_degree_of_freedom,
					Vector<Real> & lumped,
					const Vector<Int> & equation_number,
					const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_tetrahedron_10>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
void
FEMTemplate<IntegratorGauss,ShapeLagrange>::assembleLumpedTemplate<_quadrangle_8>(const Vector<Real> & field_1,
										  UInt nb_degree_of_freedom,
										  Vector<Real> & lumped,
										  const Vector<Int> & equation_number,
										  const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_quadrangle_8>(field_1, nb_degree_of_freedom,lumped, equation_number, ghost_type);
}

/* -------------------------------------------------------------------------- */
/* template instanciation                                                     */
/* -------------------------------------------------------------------------- */
//template class FEMTemplate<IntegratorGauss,ShapeLagrange>;
//template class FEMTemplate<IntegratorGauss,ShapeLinked>;

__END_AKANTU__

