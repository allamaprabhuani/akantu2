/**
 * @file   material_elastic_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Wed Aug 04 10:58:42 2010
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
inline void MaterialElastic<spatial_dimension>::computeStressOnQuad(const Matrix<Real> & grad_u,
								    Matrix<Real> & sigma,
								    const Real sigma_th) {
  Real trace = grad_u.trace();/// trace = (\nabla u)_{kk}

  /// \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla u_{ij} + \nabla u_{ji})
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      sigma(i, j) =  (i == j)*lambda*trace + mu*(grad_u(i, j) + grad_u(j, i)) + (i == j) * sigma_th;
    }
  }

}

/* -------------------------------------------------------------------------- */
template<>
inline void MaterialElastic<1>::computeStressOnQuad(const Matrix<Real> & grad_u,
						    Matrix<Real> & sigma,
						    Real sigma_th) {
  sigma(0, 0) = this->E * grad_u(0, 0) + sigma_th;
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialElastic<spatial_dimension>::computeTangentModuliOnQuad(Matrix<Real> & tangent) {
  UInt n = tangent.cols();

  //Real Ep = E/((1+nu)*(1-2*nu));
  Real Miiii = lambda + 2*mu;
  Real Miijj = lambda;
  Real Mijij = mu;

  if(spatial_dimension == 1)
    tangent(0, 0) = this->E;
  else
    tangent(0, 0) = Miiii;

  // test of dimension should by optimized out by the compiler due to the template
  if(spatial_dimension >= 2) {
    tangent(1, 1) = Miiii;
    tangent(0, 1) = Miijj;
    tangent(1, 0) = Miijj;

    tangent(n - 1, n - 1) = Mijij;
  }

  if(spatial_dimension == 3) {
    tangent(2, 2) = Miiii;
    tangent(0, 2) = Miijj;
    tangent(1, 2) = Miijj;
    tangent(2, 0) = Miijj;
    tangent(2, 1) = Miijj;

    tangent(3, 3) = Mijij;
    tangent(4, 4) = Mijij;
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void MaterialElastic<1>::computeTangentModuliOnQuad(Matrix<Real> & tangent) {
  tangent(0, 0) = E;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialElastic<spatial_dimension>::getStableTimeStep(Real h,
								  __attribute__ ((unused)) const Element & element) {
  return (h/getPushWaveSpeed());
}
