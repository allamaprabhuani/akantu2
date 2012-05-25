/**
 * @file   material_vreepeerlings_inline_impl.cc
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 * @date   Fri Feb 17 14:00:00 2012
 *
 * @brief  Specialization of the material class for the VreePeerlings material
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
inline void MaterialVreePeerlings<spatial_dimension>::computeStress(Real * F,
								    Real * sigma,
								    Real & dam,
								    Real & Equistrain,
								    Real & Kapaq) {

  types::Matrix Ft(F, 3, 3);

  Real I1=0.;
  Real J2=0.;

  if(this->plane_stress) {
     I1 = (Ft(0,0) + Ft(1,1))*(1 - 2*this->nu)/(1 - this->nu);
     Real tmp = this->nu/(this->nu - 1);
     tmp *= tmp;
     J2 = .5*(Ft(0,0)*Ft(0,0) +
	      Ft(1,1)*Ft(1,1) +
	      tmp*(Ft(0,0) + Ft(1,1))*(Ft(0,0) + Ft(1,1)) +
	      .5*(Ft(0,1) + Ft(1,0))*(Ft(0,1) + Ft(1,0))) -
       I1*I1/6.;
  }
  else {
     I1 = Ft.trace();
     J2 = .5*(Ft(0,0)*Ft(0,0) +
	      Ft(1,1)*Ft(1,1) +
	      Ft(2,2)*Ft(2,2) +
	      .5*(Ft(0,1) + Ft(1,0))*(Ft(0,1) + Ft(1,0)) +
	      .5*(Ft(1,2) + Ft(2,1))*(Ft(1,2) + Ft(2,1)) +
	      .5*(Ft(2,0) + Ft(0,2))*(Ft(2,0) + Ft(0,2))) -
       I1 * I1/6.;
  }

  Real tmp = (Kct - 1)*I1/(1 - 2*this->nu);
  Real nu1 = (1 + this->nu);

  Equistrain = tmp / (2*Kct) + 1./(2*Kct) * sqrt(tmp * tmp + (12 * Kct * J2)/(nu1 * nu1));

  MaterialElastic<spatial_dimension>::computeStress(F, sigma);

  if(!this->is_non_local) {
    computeDamageAndStress(sigma, dam, Equistrain, Kapaq);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialVreePeerlings<spatial_dimension>::computeDamageAndStress(Real * sigma,
									     Real & dam,
									     Real & Equistrain,
									     Real & Kapaq) {
  Real Fd = Equistrain - Kapaq;

  if (Fd > 0) {
    Kapaq = std::max(Equistrain,Kapaq);
    //dam = 1. - Kapa0/Kapaq*((1 - Alpha) + Alpha*exp(-Beta*(Kapaq - Kapa0)));
    dam = 1. - Kapa0/Kapaq * ((Alpha-Kapaq)/(Alpha - Kapa0));
    dam = std::min(dam,1.);
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

/* -------------------------------------------------------------------------- */
