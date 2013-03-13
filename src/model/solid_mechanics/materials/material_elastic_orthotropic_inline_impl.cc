/**
 * @file   material_elastic_orthotropic_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Tue May 08 13:01:18 2012
 *
 * @brief  Implementation of the inline functions of the orthotropic
 * elastic material
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
template<UInt spatial_dimension>
inline void
MaterialElasticOrthotropic<spatial_dimension>::computeStressOnQuad(Matrix<Real> & grad_u,
								   Matrix<Real> & sigma) {

  sigma.clear();
  /// \mathbf{\sigma} = \mathbf{S} \mathbf{F}
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      sigma(i,i) += (*S)(i,j) * grad_u(i,i);

  if(spatial_dimension == 2)
    sigma(0,1) = sigma(1,0) = (*S)(2,2) * (grad_u(0,1) + grad_u(1,0)) / 2;

  if(spatial_dimension == 3) {
    sigma(2,1) = sigma(1,2) = (*S)(3,3) * (grad_u(2,1) + grad_u(1,2)) / 2;
    sigma(0,2) = sigma(2,0) = (*S)(4,4) * (grad_u(0,2) + grad_u(2,0)) / 2;
    sigma(0,1) = sigma(1,0) = (*S)(5,5) * (grad_u(0,1) + grad_u(1,0)) / 2;
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialElasticOrthotropic<spatial_dimension>::getStableTimeStep(Real h, 
									     __attribute__ ((unused)) const Element & element) {
  return (h/getPushWaveSpeed());
}
