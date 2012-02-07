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


/* -------------------------------------------------------------------------- */
inline void MaterialMarigo::computeStress(Real * F, Real * sigma, Real & dam, Real & Y, Real &Ydq) {
  MaterialElastic::computeStress(F, sigma);

  Y = sigma[0]*F[0] +
    sigma[1]*F[1] +
    sigma[2]*F[2] +
    sigma[3]*F[3] +
    sigma[4]*F[4] +
    sigma[5]*F[5] +
    sigma[6]*F[6] +
    sigma[7]*F[7] +
    sigma[8]*F[8];

  Y *= 0.5 * (1 - dam);

  if(!is_non_local) {
    computeDamageAndStress(sigma, dam, Y, Ydq);
  }
}

/* -------------------------------------------------------------------------- */
inline void MaterialMarigo::computeDamageAndStress(Real * sigma, Real & dam, Real & Y, Real &Ydq) {
  Real Fd = Y - Ydq - Sd*dam;

  if (Fd > 0) dam = (Y - Ydq) / Sd;
  dam = std::min(dam,1.);

  sigma[0] *= 1-dam;
  sigma[4] *= 1-dam;
  sigma[8] *= 1-dam;
  sigma[1] *= 1-dam;
  sigma[3] *= 1-dam;
  sigma[2] *= 1-dam;
  sigma[6] *= 1-dam;
  sigma[5] *= 1-dam;
  sigma[7] *= 1-dam;
}

/* -------------------------------------------------------------------------- */
inline UInt MaterialMarigo::getNbDataToPack(const Element & element,
					    SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if(tag == _gst_smm_init_mat)
    size += sizeof(Real);

  size += MaterialDamage::getNbDataToPack(element, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline UInt MaterialMarigo::getNbDataToUnpack(const Element & element,
					      SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if(tag == _gst_smm_init_mat)
    size += sizeof(Real);

  size += MaterialDamage::getNbDataToPack(element, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void MaterialMarigo::packData(CommunicationBuffer & buffer,
				     const Element & element,
				     SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  if(tag == _gst_smm_init_mat)
    buffer << Yd_rand(element.type, _not_ghost)(element.element);

  MaterialDamage::packData(buffer, element, tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void MaterialMarigo::unpackData(CommunicationBuffer & buffer,
				       const Element & element,
				       SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if(tag == _gst_smm_init_mat)
    buffer >> Yd_rand(element.type, _ghost)(element.element);

  MaterialDamage::packData(buffer, element, tag);

  AKANTU_DEBUG_OUT();
}
