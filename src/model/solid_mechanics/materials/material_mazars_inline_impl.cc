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
  Real Fdiag[3] = { 0 };
  Real Fdiagp[3] = { 0 };

  Math::matrix33_eigenvalues(grad_u.storage(), Fdiag);

  Fdiagp[0] = std::max(0., Fdiag[0]);
  Fdiagp[1] = std::max(0., Fdiag[1]);
  Fdiagp[2] = std::max(0., Fdiag[2]);

  Ehat = sqrt(Fdiagp[0]*Fdiagp[0] + Fdiagp[1]*Fdiagp[1] + Fdiagp[2]*Fdiagp[2]);

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

#ifdef AKANTU_MAZARS_NON_LOCAL_AVERAGE_DAMAGE
  Real Fs = Ehat - K0;
  if (Fs > 0.) {
    Real damt;
    Real damc;
    damt =  1 - K0*(1 - At)/Ehat - At*(exp(-Bt*(Ehat - K0)));
    damc =  1 - K0*(1 - Ac)/Ehat - Ac*(exp(-Bc*(Ehat - K0)));

    Real Cdiag;
    Cdiag = E*(1-nu)/((1+nu)*(1-2*nu));

    Real SigDiag[3];
    SigDiag[0] = Cdiag*Fdiag[0] + lambda*(Fdiag[1] + Fdiag[2]);
    SigDiag[1] = Cdiag*Fdiag[1] + lambda*(Fdiag[0] + Fdiag[2]);
    SigDiag[2] = Cdiag*Fdiag[2] + lambda*(Fdiag[1] + Fdiag[0]);

    Real SigDiagT[3];
    for (UInt i = 0; i < 3; i++) {
      if(SigDiag[i] >= 0.) {
	SigDiagT[i] = SigDiag[i];
      } else {
	SigDiagT[i] = 0.;
      }
    }

    Real TraSigT;
    TraSigT = SigDiagT[0] + SigDiagT[1] + SigDiagT[2];

    Real FDiagT[3];
    for (UInt i = 0; i < 3; i++){
      FDiagT [i]=  (SigDiagT[i]*(1 + nu) - TraSigT*nu)/E;
    }

    Real alphat;
    alphat = (FDiagT[0]*Fdiagp[0] + FDiagT[1]*Fdiagp[1] + FDiagT[2]*Fdiagp[2])/(Ehat*Ehat);
    alphat = std::min(alphat, 1.);

    Real alphac;
    alphac = 1 - alphat;

    Real damtemp;
    damtemp = pow(alphat,beta)*damt + pow(alphac,beta)*damc;

    dam = std::max(damtemp, dam);
  }

  dam = std::min(dam,1.);
#endif // AKANTU_MAZARS_NON_LOCAL_AVERAGE_DAMAGE

  if(!this->is_non_local) {
    computeDamageAndStressOnQuad(grad_u, sigma, dam, Ehat);
  }
}

#ifdef AKANTU_MAZARS_NON_LOCAL_AVERAGE_DAMAGE
template<UInt spatial_dimension>
inline void
MaterialMazars<spatial_dimension>::computeDamageAndStressOnQuad(__attribute__((unused)) types::Matrix & grad_u,
								types::Matrix & sigma,
								Real & dam,
								__attribute__((unused)) Real & Ehat) {
#else
  template<UInt spatial_dimension>
inline void
MaterialMazars<spatial_dimension>::computeDamageAndStressOnQuad(types::Matrix & grad_u,
								types::Matrix & sigma,
								Real & dam,
								Real & Ehat) {
  Real Fdiag[3];
  Real Fdiagp[3];

  Math::matrix33_eigenvalues(grad_u.storage(), Fdiag);

  Fdiagp[0] = std::max(0., Fdiag[0]);
  Fdiagp[1] = std::max(0., Fdiag[1]);
  Fdiagp[2] = std::max(0., Fdiag[2]);
  Real Fs = Ehat - K0;
  if (Fs > 0.) {
    Real damt;
    Real damc;
    damt =  1 - K0*(1 - At)/Ehat - At*(exp(-Bt*(Ehat - K0)));
    damc =  1 - K0*(1 - Ac)/Ehat - Ac*(exp(-Bc*(Ehat - K0)));

    Real Cdiag;
    Cdiag = this->E*(1 - this->nu)/((1 + this->nu)*(1-2*this->nu));

    Real SigDiag[3];
    SigDiag[0] = Cdiag*Fdiag[0] + this->lambda*(Fdiag[1] + Fdiag[2]);
    SigDiag[1] = Cdiag*Fdiag[1] + this->lambda*(Fdiag[0] + Fdiag[2]);
    SigDiag[2] = Cdiag*Fdiag[2] + this->lambda*(Fdiag[1] + Fdiag[0]);

    Real SigDiagT[3];
    for (UInt i = 0; i < 3; i++) {
      if(SigDiag[i] >= 0.) {
	SigDiagT[i] = SigDiag[i];
      } else {
	SigDiagT[i] = 0.;
      }
    }

    Real TraSigT;
    TraSigT = SigDiagT[0] + SigDiagT[1] + SigDiagT[2];

    Real FDiagT[3];
    for (UInt i = 0; i < 3; i++){
      FDiagT [i]=  (SigDiagT[i]*(1 + this->nu) - TraSigT*this->nu)/this->E;
    }

    Real alphat;
    alphat = (FDiagT[0]*Fdiagp[0] + FDiagT[1]*Fdiagp[1] + FDiagT[2]*Fdiagp[2])/(Ehat*Ehat);
    alphat = std::min(alphat, 1.);

    Real alphac;
    alphac = 1 - alphat;

    Real damtemp;
    damtemp = pow(alphat,beta)*damt + pow(alphac,beta)*damc;

    dam = std::max(damtemp, dam);
  }

  dam = std::min(dam,1.);
#endif // AKANTU_MAZARS_NON_LOCAL_AVERAGE_DAMAGE

  sigma *= 1 - dam;
}
