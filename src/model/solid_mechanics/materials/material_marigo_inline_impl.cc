/**
 * @file   material_marigo_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material marigo
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
template<UInt spatial_dimension>
inline void
MaterialMarigo<spatial_dimension>::computeDamageAndStressOnQuad(types::Matrix & sigma,
								Real & dam,
								Real & Y,
								Real &Ydq) {
  Real Fd = Y - Ydq - Sd * dam;

  if (Fd > 0) dam = (Y - Ydq) / Sd;
  dam = std::min(dam,1.);

  sigma *= 1-dam;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialMarigo<spatial_dimension>::computeStressOnQuad(types::Matrix & grad_u,
						       types::Matrix & sigma,
						       Real & dam,
						       Real & Y,
						       Real &Ydq) {
  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  Y = 0;
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      Y += sigma(i,j) * grad_u(i,j);
    }
  }
  Y *= 0.5;

  //Y *= (1 - dam);
  //Y = std::min(Y, Yc);

  if(!this->is_non_local) {
    computeDamageAndStressOnQuad(sigma, dam, Y, Ydq);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline UInt MaterialMarigo<spatial_dimension>::getNbDataToPack(const Element & element,
							       SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if(tag == _gst_smm_init_mat)
    size += sizeof(Real);

  size += MaterialDamage<spatial_dimension>::getNbDataToPack(element, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline UInt MaterialMarigo<spatial_dimension>::getNbDataToUnpack(const Element & element,
								 SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if(tag == _gst_smm_init_mat)
    size += sizeof(Real);

  size += MaterialDamage<spatial_dimension>::getNbDataToPack(element, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialMarigo<spatial_dimension>::packData(CommunicationBuffer & buffer,
							const Element & element,
							SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  if(tag == _gst_smm_init_mat)
    buffer << Yd_rand(element.type, _not_ghost)(element.element);

  MaterialDamage<spatial_dimension>::packData(buffer, element, tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialMarigo<spatial_dimension>::unpackData(CommunicationBuffer & buffer,
							  const Element & element,
							  SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if(tag == _gst_smm_init_mat)
    buffer >> Yd_rand(element.type, _ghost)(element.element);

  MaterialDamage<spatial_dimension>::packData(buffer, element, tag);

  AKANTU_DEBUG_OUT();
}
