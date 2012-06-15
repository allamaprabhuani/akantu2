/**
 * @file   material_mazars_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material damage
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
MaterialMazars<spatial_dimension>::computeStressOnQuad(types::Matrix & grad_u,
						       types::Matrix & sigma,
						       Real & dam,
						       Real & Ehat) {
  types::RVector Fdiag(3);
  Fdiag.clear();
  types::RVector Fdiagp(3);
  Fdiagp.clear();

  Math::eigenvalues<spatial_dimension>(grad_u.storage(), Fdiag.storage());

  Fdiagp(0) = std::max(0., Fdiag(0));
  Fdiagp(1) = std::max(0., Fdiag(1));
  Fdiagp(2) = std::max(0., Fdiag(2));

  Ehat = sqrt(Fdiagp(0)*Fdiagp(0) + Fdiagp(1)*Fdiagp(1) + Fdiagp(2)*Fdiagp(2));

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  if(damage_in_compute_stress) {
    computeDamageOnQuad(Ehat, Fdiag, dam);
  }

  if(!this->is_non_local) {
    computeDamageAndStressOnQuad(grad_u, sigma, dam, Ehat);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialMazars<spatial_dimension>::computeDamageAndStressOnQuad(types::Matrix & grad_u,
								types::Matrix & sigma,
								Real & dam,
								Real & Ehat) {
  if(!damage_in_compute_stress) {
    types::RVector Fdiag(3);
    Fdiag.clear();
    Math::eigenvalues<spatial_dimension>(grad_u.storage(), Fdiag.storage());

    computeDamageOnQuad(Ehat, Fdiag, dam);
  }

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialMazars<spatial_dimension>::computeDamageOnQuad(const Real & Ehat,
						       const types::RVector & Fdiag,
						       Real & dam) {
  Real Fs = Ehat - K0;
  if (Fs > 0.) {
    types::RVector Fdiagp(3);
    Fdiagp.clear();

    Fdiagp(0) = std::max(0., Fdiag(0));
    Fdiagp(1) = std::max(0., Fdiag(1));
    Fdiagp(2) = std::max(0., Fdiag(2));

    Real damt;
    Real damc;
    damt = 1 - K0*(1 - At)/Ehat - At*(exp(-Bt*(Ehat - K0)));
    damc = 1 - K0*(1 - Ac)/Ehat - Ac*(exp(-Bc*(Ehat - K0)));

    Real Cdiag;
    Cdiag = this->E*(1-this->nu)/((1+this->nu)*(1-2*this->nu));
    
    types::RVector SigDiag(3);
    SigDiag(0) = Cdiag*Fdiag(0) + this->lambda*(Fdiag(1) + Fdiag(2));
    SigDiag(1) = Cdiag*Fdiag(1) + this->lambda*(Fdiag(0) + Fdiag(2));
    SigDiag(2) = Cdiag*Fdiag(2) + this->lambda*(Fdiag(1) + Fdiag(0));

    types::RVector SigDiagT(3);
    for (UInt i = 0; i < 3; i++) {
      if(SigDiag(i) >= 0.) {
	SigDiagT(i) = SigDiag(i);
      } else {
	SigDiagT(i) = 0.;
      }
    }

    Real TraSigT;
    TraSigT = SigDiagT(0) + SigDiagT(1) + SigDiagT(2);

    types::RVector FDiagT(3);
    for (UInt i = 0; i < 3; i++){
      FDiagT (i)=  (SigDiagT(i)*(1 + this->nu) - TraSigT*this->nu)/this->E;
    }

    Real alphat;
    alphat = (FDiagT(0)*Fdiagp(0) + FDiagT(1)*Fdiagp(1) + FDiagT(2)*Fdiagp(2))/(Ehat*Ehat);
    alphat = std::min(alphat, 1.);

    Real alphac;
    alphac = 1 - alphat;

    Real damtemp;
    damtemp = pow(alphat,beta)*damt + pow(alphac,beta)*damc;

    dam = std::max(damtemp, dam);
  }

  dam = std::min(dam,1.);
}
