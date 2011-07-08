/**
 * @file   material_neohookean_inline_impl.cc
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material neohookean
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
#include "aka_math.hh"

/* -------------------------------------------------------------------------- */
inline void MaterialNeohookean::computeStress(Real * F, Real * sigma) {
  ///First compute the Left Cauchy-Green deformation tensor : C= F^tF.
  Real C[3*3];
  Math::matMul<true,false>(3,3,3,1.,F,F,0.,C);

  ///Compute determinant of C
  Real detC;
  detC=Math::det3(C);

  Real defvol ;
  defvol= 0.5*log(detC);

  Real p;
  p = lambda*defvol;

  Real traceC;
  traceC = C[0]+C[4]+C[8];

  Real Cinv[3*3];
  Math::inv3(C, Cinv);

  Real coef = p - mu;

  Real S[3*3];
  for(UInt i=0; i < 3*3; i++){
    S[i] = coef*Cinv[i];
  }

  S[0] = S[0] + mu;
  S[4] = S[4] + mu;
  S[8] = S[8] + mu;

  Math::matrix_matrix(3, 3,3, F, S, sigma, 1.);
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void  MaterialNeohookean::computeTangentStiffness(Real * tangent, Real * F ) {
  UInt n = (dim * (dim - 1) / 2 + dim);
  Real J = Math::det3(F);
  Real alpha = mu - 0.5*lambda*(J*J - 1);
  Real Miiii = 2*mu+lambda;
  Real Miijj = lambda*J*J;
  Real Mijij = alpha;

  tangent[0 * n + 0] = Miiii;

  // test of dimension should by optimized out by the compiler due to the template
  if(dim >= 2) {
    tangent[1 * n + 1] = Miiii;
    tangent[0 * n + 1] = Miijj;
    tangent[1 * n + 0] = Miijj;

    tangent[(n - 1) * n + (n - 1)] = Mijij;
  }

  if(dim == 3) {
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
inline void MaterialNeohookean::computePotentialEnergy(Real * F, Real * epot) {
  *epot = 0.;
  Real C[3*3];
  Math::matMul<true,false>(3,3,3,1.,F,F,0.,C);

  ///Compute determinant of C
  Real detC;
  detC = Math::det3(C);

  Real defvol ;
  defvol = 0.5*log(detC);

  Real p;
  p = lambda*defvol;

  Real traceC;
  traceC = C[0] + C[4] + C[8];

  ///energie potentielle
  *epot = (0.5*p - mu)*defvol + 0.5*mu*(traceC - 3.);
}

/* -------------------------------------------------------------------------- */
inline Real MaterialNeohookean::getStableTimeStep(Real h,
						 const Element & element) {
  return (h/celerity(element));
}
