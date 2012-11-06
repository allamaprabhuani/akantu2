/**
 * @file   material_vreepeerlings_inline_impl.cc
 *
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date   Fri Feb 24 14:27:15 2012
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
inline void
MaterialVreePeerlings<spatial_dimension>::computeStressOnQuad(types::RMatrix & grad_u,
							      types::RMatrix & sigma,
							      Real & dam,
							      Real & Equistrain,
							      Real & Equistrain_rate,
							      Real & Kapaq,
							      Real dt,
							      types::RMatrix & strain_rate_vrplgs,
							      Real & crit_strain) {
  Real I1=0.;
  Real J2=0.;
  Real I1rate=0.;
  Real J2rate=0.;

  if(this->plane_stress) {
     I1 = (grad_u(0,0) + grad_u(1,1))*(1 - 2*this->nu)/(1 - this->nu);
     I1rate = (strain_rate_vrplgs(0,0) + strain_rate_vrplgs(1,1))*(1 - 2*this->nu)/(1 - this->nu);

     Real tmp = this->nu/(this->nu - 1);
     tmp *= tmp;

     J2 = .5*(grad_u(0,0)*grad_u(0,0) +
	      grad_u(1,1)*grad_u(1,1) +
	      tmp*(grad_u(0,0) + grad_u(1,1))*(grad_u(0,0) + grad_u(1,1)) +
	      .5*(grad_u(0,1) + grad_u(1,0))*(grad_u(0,1) + grad_u(1,0))) -
       I1*I1/6.;
     //     J2rate = .5*(strain_rate_vrplgs(0,0)*strain_rate_vrplgs(0,0) +
     //  	      strain_rate_vrplgs(1,1)*strain_rate_vrplgs(1,1) +
     //               tmp*(strain_rate_vrplgs(0,0) + strain_rate_vrplgs(1,1))*(strain_rate_vrplgs(0,0) + strain_rate_vrplgs(1,1)) +
     //	             .5*(strain_rate_vrplgs(0,1) + strain_rate_vrplgs(1,0))*(strain_rate_vrplgs(0,1) + strain_rate_vrplgs(1,0))) -
     //  I1rate*I1rate/6.;

     J2rate = (grad_u(0,0)*strain_rate_vrplgs(0,0) +
  	      grad_u(1,1)*strain_rate_vrplgs(1,1) +
              tmp*(grad_u(0,0) + grad_u(1,1))*(strain_rate_vrplgs(0,0) + strain_rate_vrplgs(1,1)) +
	      (grad_u(0,1)*strain_rate_vrplgs(1,0)) + (grad_u(0,1)*strain_rate_vrplgs(1,0))) -
       I1*I1rate/3.;

  }
  else {
     I1 = grad_u.trace();
     I1rate = strain_rate_vrplgs.trace();

     J2 = .5*(grad_u(0,0)*grad_u(0,0) +
	      grad_u(1,1)*grad_u(1,1) +
	      grad_u(2,2)*grad_u(2,2) +
	      .5*(grad_u(0,1) + grad_u(1,0))*(grad_u(0,1) + grad_u(1,0)) +
	      .5*(grad_u(1,2) + grad_u(2,1))*(grad_u(1,2) + grad_u(2,1)) +
	      .5*(grad_u(2,0) + grad_u(0,2))*(grad_u(2,0) + grad_u(0,2))) -
       I1 * I1/6.;

     // J2rate = .5*(strain_rate_vrplgs(0,0)*strain_rate_vrplgs(0,0) +
     //	      strain_rate_vrplgs(1,1)*strain_rate_vrplgs(1,1) +
     //	      strain_rate_vrplgs(2,2)*strain_rate_vrplgs(2,2) +
     //	      .5*(strain_rate_vrplgs(0,1) + strain_rate_vrplgs(1,0))*(strain_rate_vrplgs(0,1) + strain_rate_vrplgs(1,0)) +
     //	      .5*(strain_rate_vrplgs(1,2) + strain_rate_vrplgs(2,1))*(strain_rate_vrplgs(1,2) + strain_rate_vrplgs(2,1)) +
     //	      .5*(strain_rate_vrplgs(2,0) + strain_rate_vrplgs(0,2))*(strain_rate_vrplgs(2,0) + strain_rate_vrplgs(0,2))) -
     //  I1rate * I1rate/6.;

     J2rate = (grad_u(0,0)*strain_rate_vrplgs(0,0) +
    	      grad_u(1,1)*strain_rate_vrplgs(1,1) +
    	      grad_u(2,2)*strain_rate_vrplgs(2,2) +
    	      (grad_u(0,1)*strain_rate_vrplgs(1,0)) + (grad_u(0,1)*strain_rate_vrplgs(1,0)) +
    	      (grad_u(1,2)*strain_rate_vrplgs(2,1)) + (grad_u(1,2)*strain_rate_vrplgs(2,1)) +
    	      (grad_u(2,0)*strain_rate_vrplgs(0,2)) + (grad_u(2,0)*strain_rate_vrplgs(0,2))) -
      I1 * I1rate/6.;
  }

  //Real tmp = (Kct - 1)*I1/(1 - 2*this->nu);
  Real tmp = (Kct - 1)/(1 - 2*this->nu);
  Real tmp2 = (12 * Kct)/((1 + this->nu)*(1 + this->nu));
  //  Real tmp3 = (12)/((1 + this->nu)*(1 + this->nu));
  ///////////////////////////////////////////////////debut
  //Real tmprate = (Kct - 1)*I1rate/(1 - 2*this->nu);
  ///////////////////////////////////////////////////fin
  // Real nu1 = (1 + this->nu);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////debut
  
  Real Equistrainratenew=tmp*I1rate / (2*Kct) + 1./(4*Kct) * (2*tmp*tmp*I1*I1rate + tmp2*J2rate) / sqrt(tmp*tmp*I1*I1 + tmp2*J2);
  Equistrain = tmp * I1 / (2*Kct) + 1./(2*Kct) * sqrt(tmp*tmp*I1*I1 + tmp2*J2);
 
  if(Equistrainratenew*Equistrain>0){
    if(I1rate * std::abs(I1rate) > 0){
      Equistrain_rate = Equistrainratenew;
    }else{
      Equistrain_rate =tmp * std::abs(I1rate) / (2*Kct) + 1./(4*Kct)*(2*tmp*tmp*I1*I1rate + tmp2*J2rate) / sqrt(tmp*tmp*I1*I1 + tmp2*J2);
    }
    //test
    //Equistrain_rate = 1./(4.) * (tmp3*J2rate) / sqrt(tmp3*J2);
  }

  // Real tmprate = (Kct - 1)*abs(I1rate)/(1 - 2*this->nu);
  // Real Equistrainratenew = tmprate / (2*Kct) + 1./(2*Kct) * sqrt(tmprate * tmprate + (12 * Kct * J2rate)/(nu1 * nu1));
  //   if((I1*abs(I1)>0. && I1rate*abs(I1rate)>0.) || (I1*abs(I1)<0. && I1rate*abs(I1rate)<0.)){
  //    Equistrain_rate=Equistrainratenew;
  //  }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////fin

  // Equistrain = tmp / (2*Kct) + 1./(2*Kct) * sqrt(tmp * tmp + (12 * Kct * J2)/(nu1 * nu1));
  // Equistrain = tmp * I1 / (2*Kct) + 1./(2*Kct) * sqrt(tmp*tmp * I1*I1 + tmp2 * J2);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////debut
     // Equistrain_rate = tmprate / (2*Kct) + 1./(2*Kct) * sqrt(tmprate * tmprate + (12 * Kct * J2rate)/(nu1 * nu1));
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////fin

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  if(!this->is_non_local) {
    computeDamageAndStressOnQuad(sigma, dam, Equistrain, Equistrain_rate, Kapaq, dt, crit_strain);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialVreePeerlings<spatial_dimension>::computeDamageAndStressOnQuad(types::RMatrix & sigma,
								       Real & dam,
								       Real & Equistrain,
								       Real & Equistrain_rate,
								       Real & Kapaq,
								       __attribute__((unused)) Real dt,
								       Real & crit_strain) {
  ///Richeton Formula
  // Real cst1 = 190e6 - 470e3 * 298.;
  // Real cst2 = 2*1.38*1000000*298/5.14;
  // Real cst3 = 7.46e15*exp(-1*((90000)/(6.022*1.38*298)));
  // Real cst4 = Equistrain_rate/cst3;
  // Real cst5 = pow(cst4,1/6.37);
  // Real cst6 = log(cst5 + sqrt(1+cst5*cst5));
  // Real cst7 = cst1+cst2*cst6;
  // Real cst8 = cst7/(3090.*1e6*1.48);
  ///Zhou & Molinari Formula 
  //  Real cst10 = (75./3090.)*(1.+(Equistrain_rate/Beta)); 
  // Real cst10 = (Kapa0)*(1.+(Equistrain_rate/Beta)); 
  /// Real cst10 = (75./3090.)*(1.+(Equistrain_rate/(pow(Beta,0.98+1./(1.+Equistrain_rate))))); 
  // Real cst11 = Alpha;
  // Real cst11 = Alpha*(1.+(Equistrain_rate/Beta)); 
  //  Real cst11 = 0;
  //  if (Equistrain_rate >= 750000.){
  // cst11 = Alpha*(1.+(Equistrain_rate/500000.)); 
  //} else {
  // cst11 = Alpha*(1.+(Equistrain_rate/1000000.)); 
  //}
 
  // epsilono
  //   Real cst10 = Kapa0;
  //     Real cst10 = (Kapa0)*(1.+(Equistrain_rate/Beta)); 
       Real cst10 = (75./3090.)*(1.+(Equistrain_rate/Beta)); 
  //    Real cst10 = Equistrain_rate; 

  //test
  //       Real sigmaco=(this->E)*Kapa0*Alpha/(Alpha-Kapa0);
  //    Real cst10= (sigmaco/this->E)*(1.+(Equistrain_rate/Beta))/(1.+sigmaco/(Alpha*this->E));
  //deg2
  //    Real cst10 = (Kapa0)*(1.+(Equistrain_rate/Beta)+(Equistrain_rate/Beta)*(Equistrain_rate/Beta)); 
  //ln
  // Real cst10 = (Kapa0)*(1.+(log(Equistrain_rate+1.)/(Beta/10000.))); 

  // epsilonc
         Real cst11 = Alpha;
  //  Real cst11 = Alpha*(1.+(Equistrain_rate/Beta)); 



  Real cst9 = 0;
  if (cst10 > Kapa0) {
    if (cst10 < 0.16) {
      cst9=cst10;
    } else {
      cst9=0.16;
    }
  } else {
    cst9 = Kapa0;
  }

    Real Fd1 = Equistrain - cst9; 
    if (Fd1 > 0)
     {
       //Real Fd = Equistrain - Kapaq;
       // if (Fd  > 0)
       //{
	      Real dam1 = 1. - cst9/Equistrain * ((cst11-Equistrain)/(cst11 - cst9));
	   //  Real dam1 = 1. - cst9/Equistrain * ((Alpha-Equistrain)/(Alpha - cst9));

        if (dam1 > dam){
           Kapaq = std::max(Equistrain,Kapaq);
           if (dam1 >= 0.9999){
  	     dam =1.;
  	   }
  	   else {    
  	     dam = std::min(dam1,1.);
             crit_strain=cst9;
  	   }
	 }
        //}
     }

  // Real Alpharates= Alpha/(1.+Equistrain_rate/(6000000.)); new law 2
  //   Real Alpharates= 0.6/(1.+Equistrain_rate/(1000000.));
  //new law 3
  // Real Alpharates = 0.6/(1+Equistrain_rate/(pow(1e6,0.98+1./(1+Equistrain_rate))));
  //  Real Alpharate=0;
  //    if (Alpharates <= 0.15){
  //      Alpharate=0.15;
  //    }
  //    else {
  //     Alpharate=Alpharates;
  //    }
  //  Real dam2 = 1. - cst9/Equistrain * ((Alpharate - Equistrain)/(Alpharate - cst9));
  //  if (dam2 >= dam){
  //  dam = std::min(dam2,1.);
  //  }
   sigma *= 1-dam;
   //   crit_strain=cst9;
}

/* -------------------------------------------------------------------------- */
