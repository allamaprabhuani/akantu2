/**
 * @file   fem_template_cohesive_tmpl.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Nov  5 17:23:44 2012
 *
 * @brief  Cohesive element template specialization declaration
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


/* -------------------------------------------------------------------------- */
template <>
void FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> >::
gradientOnQuadraturePoints(const Array<Real> &u,
			   Array<Real> &nablauq,
			   const UInt nb_degree_of_freedom,
			   const ElementType & type,
			   const GhostType & ghost_type,
			   const Array<UInt> * filter_elements) const;

/* -------------------------------------------------------------------------- */
template <>
void FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> >::
initShapeFunctions(const GhostType & ghost_type);

/* -------------------------------------------------------------------------- */
template <>
Real FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> >::
integrate(const Array<Real> & f,
	  const ElementType & type,
	  const GhostType & ghost_type,
	  const Array<UInt> * filter_elements) const;
/* -------------------------------------------------------------------------- */
template <>
void FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> >::
integrate(const Array<Real> & f,
	  Array<Real> &intf,
	  UInt nb_degree_of_freedom,
	  const ElementType & type,
	  const GhostType & ghost_type,
	  const Array<UInt> * filter_elements) const;
