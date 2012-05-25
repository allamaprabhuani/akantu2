/**
 * @file   material_elastic_orthotropic_inline_impl.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Apr 12 13:40:42 2012
 *
 * @brief Implementation of the inline functions of the orthotropic
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
inline void MaterialElasticOrthotropic<spatial_dimension>::computeStress(Real * F, Real * sigma) {
  /// \mathbf{\sigma} = \mathbf{S} \mathbf{F}
  sigma[0] = S[0] * F[0] + S[3] * F[4] + S[4] * F[8];
  sigma[4] = S[3] * F[0] + S[1] * F[4] + S[5] * F[8];
  sigma[8] = S[4] * F[0] + S[5] * F[4] + S[2] * F[8];

  sigma[1] = sigma[3] = S[6] * (F[1] + F[3]) / 2;
  sigma[2] = sigma[6] = S[7] * (F[2] + F[6]) / 2;
  sigma[5] = sigma[7] = S[8] * (F[5] + F[7]) / 2;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialElasticOrthotropic<dim>::computeTangentStiffness(Real * tangent) {

  UInt n = (dim * (dim - 1) / 2 + dim);

  tangent[0 * n + 0] = S[0];

  // test of dimension should by optimized out by the compiler due to the template
  if(dim >= 2) {
    tangent[1 * n + 1] = S[1];
    tangent[0 * n + 1] = S[3];
    tangent[1 * n + 0] = S[3];

    tangent[(n - 1) * n + (n - 1)] = S[6];
  }

  if(dim == 3) {
    tangent[2 * n + 2] = S[2];
    tangent[0 * n + 2] = S[4];
    tangent[1 * n + 2] = S[5];
    tangent[2 * n + 0] = S[4];
    tangent[2 * n + 1] = S[5];

    tangent[3 * n + 3] = S[7];
    tangent[4 * n + 4] = S[8];
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialElasticOrthotropic<spatial_dimension>::getStableTimeStep(Real h, 
									     __attribute__ ((unused)) const Element & element) {
  return (h/getPushWaveSpeed());
}
