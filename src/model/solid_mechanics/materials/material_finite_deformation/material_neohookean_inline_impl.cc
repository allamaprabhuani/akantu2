/**
 * @file   material_neohookean_inline_impl.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
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

#include <cmath>
#include "material_neohookean.hh"


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialNeohookean<dim>::computeDeltaStressOnQuad(const Matrix<Real> & grad_u, const Matrix<Real> & grad_delta_u,
							   Matrix<Real> & delta_S){


}

template<UInt dim>
inline void MaterialNeohookean<dim>::computeStressOnQuad(Matrix<Real> & grad_u,
						      Matrix<Real> & sigma) {
  //Neo hookean book
  Matrix<Real> F(dim, dim);
  Matrix<Real> C(dim, dim);//Right green
  Matrix<Real> Cminus(dim, dim);//Right green


  this->template gradUToF<dim > (grad_u, F);
  this->rightCauchy(F, C);
  Real J = F.det();
  Cminus.inverse(C);

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      sigma(i, j) = (i == j) *  mu  + (lambda * log(J) - mu) * Cminus(i, j);

}

template<UInt dim>
inline void MaterialNeohookean<dim>::computePiolaKirchhoffOnQuad(const Matrix<Real> & E,
							      Matrix<Real> & S) {

  Real trace = E.trace(); /// trace = (\nabla u)_{kk}

  /// \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla u_{ij} + \nabla u_{ji})
  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      S(i, j) = (i == j) * lambda * trace + 2.0 * mu * E(i, j);

}


/**************************************************************************************/
/*  Computation of the potential energy for a this neo hookean material */
template<UInt dim>
inline void MaterialNeohookean<dim>::computePotentialEnergyOnQuad(const Matrix<Real> & grad_u,
                                                               Real & epot){
  Matrix<Real> F(dim, dim);
  Matrix<Real> C(dim, dim);//Right green

  this->template gradUToF<dim > (grad_u, F);
  this->rightCauchy(F, C);
  Real J = F.det();

  epot=0.5*lambda*pow(log(J),2.)+ mu * (-log(J)+0.5*(C.trace()-dim));
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialNeohookean<dim>::computeTangentModuliOnQuad(Matrix<Real> & tangent, Matrix<Real> & grad_u) {

  //Neo hookean book
  UInt cols = tangent.cols();
  UInt rows = tangent.rows();
  Matrix<Real> F(dim, dim);
  Matrix<Real> C(dim, dim);
  Matrix<Real> Cminus(dim, dim);
  this->template gradUToF<dim > (grad_u, F);
  this->rightCauchy(F, C);
  Real J = F.det();
  Cminus.inverse(C);

  for (UInt m = 0; m < rows; m++) {
    UInt i = VoigtHelper<dim>::vec[m][0];
    UInt j = VoigtHelper<dim>::vec[m][1];
    for (UInt n = 0; n < cols; n++) {
      UInt k = VoigtHelper<dim>::vec[n][0];
      UInt l = VoigtHelper<dim>::vec[n][1];

      //book belytchko
      tangent(m, n) = lambda * Cminus(i, j) * Cminus(k, l) +
	(mu - lambda * log(J)) * (Cminus(i, k) * Cminus(j, l) + Cminus(i, l) * Cminus(k, j));

    }
  }

}

/* -------------------------------------------------------------------------- */
/*template<>
  inline void MaterialNeohookean < 1 > ::computeTangentModuliOnQuad(Matrix<Real> & tangent) {
  tangent(0, 0) = E;
  }*/

