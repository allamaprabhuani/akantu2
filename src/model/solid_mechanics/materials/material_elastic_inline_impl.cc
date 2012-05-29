/**
 * @file   material_elastic_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material elastic
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
inline void MaterialElastic<spatial_dimension>::computeStress(Real * F, Real * sigma) {
  Real trace = F[0] + F[4] + F[8]; /// \F_{11} + \F_{22} + \F_{33}

  /// \sigma_{ij} = \lamda * \F_{kk} * \delta_{ij} + 2 * \mu * \F_{ij}
  sigma[0] = lambda * trace + 2*mu*F[0];
  sigma[4] = lambda * trace + 2*mu*F[4];
  //  if(plane_stress) F[8] = (F[0] + F[4])*(nu/(nu-1.));
  sigma[8] = lambda * trace + 2*mu*F[8];


  sigma[1] = sigma[3] =  mu * (F[1] + F[3]);
  sigma[2] = sigma[6] =  mu * (F[2] + F[6]);
  sigma[5] = sigma[7] =  mu * (F[5] + F[7]);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialElastic<spatial_dimension>::computeStress(types::Matrix & grad_u,
							      types::Matrix & sigma) {
  Real trace = grad_u.trace();/// trace = (\nabla u)_{kk}

  /// \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla u_{ij} + \nabla u_{ji})
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      sigma(i, j) =  (i == j) * lambda * trace + mu*(grad_u(i, j) + grad_u(j, i));
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void MaterialElastic<1>::computeStress(Real * F, Real * sigma) {
  sigma[0] = E * F[0];
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::computeTangentStiffness(Real * tangent) {

  UInt n = (spatial_dimension * (spatial_dimension - 1) / 2 + spatial_dimension);

  Real Ep = E/((1+nu)*(1-2*nu));
  Real Miiii = Ep * (1-nu);
  Real Miijj = Ep * nu;
  Real Mijij = Ep * (1-2*nu) * .5;

  tangent[0 * n + 0] = Miiii;

  // test of dimension should by optimized out by the compiler due to the template
  if(spatial_dimension >= 2) {
    tangent[1 * n + 1] = Miiii;
    tangent[0 * n + 1] = Miijj;
    tangent[1 * n + 0] = Miijj;

    tangent[(n - 1) * n + (n - 1)] = Mijij;
  }

  if(spatial_dimension == 3) {
    tangent[2 * n + 2] = Miiii;
    tangent[0 * n + 2] = Miijj;
    tangent[1 * n + 2] = Miijj;
    tangent[2 * n + 0] = Miijj;
    tangent[2 * n + 1] = Miijj;

    tangent[3 * n + 3] = Mijij;
    tangent[4 * n + 4] = Mijij;
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialElastic<spatial_dimension>::getStableTimeStep(Real h,
								  __attribute__ ((unused)) const Element & element) {
  return (h/getPushWaveSpeed());
}
