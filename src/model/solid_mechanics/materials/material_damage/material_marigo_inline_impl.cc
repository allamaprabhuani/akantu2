/**
 * @file   material_marigo_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Feb 02 11:09:36 2012
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
MaterialMarigo<spatial_dimension>::computeStressOnQuad(types::RMatrix & grad_u,
						       types::RMatrix & sigma,
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

  if(damage_in_y) Y *= (1 - dam);

  if(yc_limit) Y = std::min(Y, Yc);

  if(!this->is_non_local) {
    computeDamageAndStressOnQuad(sigma, dam, Y, Ydq);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialMarigo<spatial_dimension>::computeDamageAndStressOnQuad(types::RMatrix & sigma,
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
inline UInt MaterialMarigo<spatial_dimension>::getNbDataForElements(const Vector<Element> & elements,
								    SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if(tag == _gst_smm_init_mat) {
    size += sizeof(Real) * this->getModel().getNbQuadraturePoints(elements);
  }

  size += MaterialDamage<spatial_dimension>::getNbDataForElements(elements, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialMarigo<spatial_dimension>::packElementData(CommunicationBuffer & buffer,
							       const Vector<Element> & elements,
							       SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  if(tag == _gst_smm_init_mat) {
    this->packElementDataHelper(Yd_rand, buffer, elements);
    // UInt nb_quadrature_points = this->model->getFEM().getNbQuadraturePoints(element.type);
    // Vector<Real>::const_iterator<Real> Yds = Yd_rand(element.type, _not_ghost).begin();
    // Yds += element.element * nb_quadrature_points;
    // for (UInt q = 0; q < nb_quadrature_points; ++q, ++Yds)
    //   buffer << *Yds;
  }

  MaterialDamage<spatial_dimension>::packElementData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialMarigo<spatial_dimension>::unpackElementData(CommunicationBuffer & buffer,
								 const Vector<Element> & elements,
								 SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if(tag == _gst_smm_init_mat) {
    this->unpackElementDataHelper(Yd_rand, buffer, elements);
    // UInt nb_quadrature_points = this->model->getFEM().getNbQuadraturePoints(element.type);
    // Vector<Real>::iterator<Real> Ydr = Yd_rand(element.type, _ghost).begin();
    // Ydr += element.element * nb_quadrature_points;
    // for (UInt q = 0; q < nb_quadrature_points; ++q, ++Ydr)
    //   buffer << *Ydr;
  }

  MaterialDamage<spatial_dimension>::unpackElementData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}
