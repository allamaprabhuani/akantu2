/**
 * @file   material_mazars_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
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
inline void MaterialMazars::computeStress(Real * F, Real * sigma, Real & dam) {

  Real trace = F[0] + F[4] + F[8]; 
  ///  @f$\F_{11} + \F_{22} + \F_{33} @f$
  /// @f$ \sigma_{ij} = \lamda * \F_{kk} * \delta_{ij} + 2 * \mu * \F_{ij} @f$
  Real K = 1./3.*(E/(1.-2.*nu));
  Real G = E/(2*(1+nu));
  Real Fdiag[3];
  Real Fdiagp[3];

  Math::matrix33_eigenvalues(F,Fdiag);

  Fdiagp[0] = std::max(0.,Fdiag[0]);
  Fdiagp[1] = std::max(0.,Fdiag[1]);
  Fdiagp[2] = std::max(0.,Fdiag[2]);
  
  Real Ehat=sqrt(Fdiagp[0]*Fdiagp[0]+Fdiagp[1]*Fdiagp[1]+Fdiagp[2]*Fdiagp[2]);

  sigma[0] = K * trace +2*G*(F[0]-trace/3);
  sigma[4] = K * trace +2*G*(F[4]-trace/3);
  sigma[8] = K * trace +2*G*(F[8]-trace/3);
  sigma[1] = sigma[3] =  G * (F[1] + F[3]);
  sigma[2] = sigma[6] =  G * (F[2] + F[6]);
  sigma[5] = sigma[7] =  G * (F[5] + F[7]);

 
  Real Fs = Ehat-K0;
  Real damt  ;
  Real damc  ;
  if (Fs > 0) {
    damt =  1-K0*(1-At)/Ehat-At*(exp(-Bt*(Ehat-K0)));
    damc =  1-K0*(1-Ac)/Ehat-Ac*(exp(-Bc*(Ehat-K0)));
    Real Cdiag ;
    Cdiag = E*(1-nu)/((1+nu)*(1-2*nu));
    
    Real SigDiag[3];
    SigDiag[0] = Cdiag*Fdiag[0]+ lambda*(Fdiag[1]+Fdiag[2]);
    SigDiag[1] = Cdiag*Fdiag[1]+ lambda*(Fdiag[0]+Fdiag[2]);
    SigDiag[2] = Cdiag*Fdiag[2]+ lambda*(Fdiag[1]+Fdiag[0]);
    Real SigDiagT[3];
    Real SigDiagC[3];   
    for (UInt i = 0; i<3;i++){
      if( SigDiag[i]>=0.){
	 SigDiagT[i]= SigDiag[i];
	 SigDiagC[i]=0. ;}
      else{
	 SigDiagC[i]= SigDiag[i];
	 SigDiagT[i]=0. ;}	
    }
    Real TraSigT, TraSigC;
    TraSigT =SigDiagT[0]+ SigDiagT[1]+SigDiagT[2];
    TraSigC =SigDiagC[0]+ SigDiagC[1]+SigDiagC[2];
    Real FDiagT [3];
       for (UInt i = 0; i<3;i++){
	 FDiagT [i]=  (SigDiagT[i]*(1+nu)-TraSigT*nu)/E;
       }
    Real alphat ;
    Real alphac ;
    Real damtemp;
    alphat= (FDiagT[0]*Fdiagp[0]+FDiagT[1]*Fdiagp[1]+FDiagT[2]*Fdiagp[2])/(Ehat*Ehat);
    alphat=std::min(alphat,1.);
    alphac = 1-alphat;
    damtemp= pow(alphat,beta)*damt+pow(alphac,beta)*damc;
  dam=std::max(damtemp,dam);
  }
  dam = std::min(dam,1.);

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
inline void MaterialMazars::computePotentialEnergy(Real * F, Real * sigma, Real * epot) {
  *epot = 0.;
  for (UInt i = 0, t = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j, ++t)
      *epot += sigma[t] * F[t] ;
  *epot *= .5;
}

/* -------------------------------------------------------------------------- */
inline Real MaterialMazars::celerity() {
  return sqrt(E/rho);
}

/* -------------------------------------------------------------------------- */
inline Real MaterialMazars::getStableTimeStep(Real h) {
  return (h/celerity());
}
