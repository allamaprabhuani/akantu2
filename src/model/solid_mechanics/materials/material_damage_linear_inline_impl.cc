/**
 * @file   material_damage_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material damage linear
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
inline void MaterialDamageLinear::computeStress(Real * F, Real * sigma, Real & dam, Real &K) {
  //  Real trace = F[0] + F[4] + F[8];
  // Real K = 1./3. * (E/(1. - 2.*nu));
  // Real G = E / (2*(1 + nu));
  // Real lambda = nu * E / ((1 + nu) * (1 - 2*nu));


  Real Fdiag[3];
  Real Fdiagp[3];

  Math::matrix33_eigenvalues(F, Fdiag);

  Fdiagp[0] = std::max(0., Fdiag[0]);
  Fdiagp[1] = std::max(0., Fdiag[1]);
  Fdiagp[2] = std::max(0., Fdiag[2]);

  Real Ehat=sqrt(Fdiagp[0]*Fdiagp[0]+Fdiagp[1]*Fdiagp[1]+Fdiagp[2]*Fdiagp[2]);

  MaterialElastic::computeStress(F, sigma);

  Real Fd = Ehat-K;

  if (Fd > 0) {
    dam = (Ehat - Epsmin) / (Epsmax-Epsmin)*(Ehat/Epsmax);
    dam = std::min(dam, 1.);
    K=Ehat;
  }

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

