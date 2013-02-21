/**
 * @file   material_cohesive_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Feb 23 17:58:45 2012
 *
 * @brief  MaterialCohesive inline implementation
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
template<ElementType type>
void MaterialCohesive::computeNormal(const Vector<Real> & position,
				     Vector<Real> & normal,
				     GhostType ghost_type) {
}

/* -------------------------------------------------------------------------- */
inline UInt MaterialCohesive::getNbDataForElements(const Vector<Element> & elements,
						   SynchronizationTag tag) const {

  if(tag == _gst_smm_stress) {
    return spatial_dimension * sizeof(Real) * this->getModel().getNbQuadraturePoints(elements, "CohesiveFEM");
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
inline void MaterialCohesive::packElementData(CommunicationBuffer & buffer,
					      const Vector<Element> & elements,
					      SynchronizationTag tag) const {
  if(tag == _gst_smm_stress) {
    packElementDataHelper(tractions, buffer, elements, "CohesiveFEM");
  }

}

/* -------------------------------------------------------------------------- */
inline void MaterialCohesive::unpackElementData(CommunicationBuffer & buffer,
						const Vector<Element> & elements,
						SynchronizationTag tag) {
  if(tag == _gst_smm_stress) {
    unpackElementDataHelper(tractions, buffer, elements, "CohesiveFEM");
  }

}
