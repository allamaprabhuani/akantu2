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
template<UInt spatial_dimension, template <UInt> class MatParent>
inline void
MaterialVreePeerlings<spatial_dimension, MatParent>::computeStressOnQuad(Matrix<Real> & grad_u,
                                                                         Matrix<Real> & sigma,
                                                                         Real & dam,
                                                                         Real & Equistrain,
                                                                         Real & Equistrain_rate,
                                                                         Real & Kapaq,
                                                                         Real dt,
                                                                         Matrix<Real> & strain_rate_vrplgs,
                                                                         Real & crit_strain,
                                                                         Real & crit_strain_rate,
                                                                         Real & rdr_damage,
                                                                         Real & nb_damage) {
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
  
  // Real Equistrainratenew=tmp*I1rate / (2*Kct) + 1./(4*Kct) * (2*tmp*tmp*I1*I1rate + tmp2*J2rate) / sqrt(tmp*tmp*I1*I1 + tmp2*J2);
  Equistrain = tmp * I1 / (2*Kct) + 1./(2*Kct) * sqrt(tmp*tmp*I1*I1 + tmp2*J2);
 
  //  if(Equistrainratenew*Equistrain>0){
    if(I1 * I1rate > 0 || J2rate > 0){
      Equistrain_rate =tmp * std::abs(I1rate) / (2*Kct) + 1./(4*Kct)*(2*tmp*tmp*I1*I1rate + tmp2*J2rate) / sqrt(tmp*tmp*I1*I1 + tmp2*J2);
    }else{
      Equistrain_rate = Equistrain_rate;
    }
    //test
    //Equistrain_rate = 1./(4.) * (tmp3*J2rate) / sqrt(tmp3*J2);
    //  }

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

  //  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  if(!this->is_non_local) {
    computeDamageAndStressOnQuad(sigma, dam, Equistrain, Equistrain_rate, Kapaq, dt, crit_strain, crit_strain_rate, rdr_damage, nb_damage);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class MatParent>
inline void
MaterialVreePeerlings<spatial_dimension, MatParent>::computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                                                                  Real & dam,
                                                                                  Real & Equistrain,
                                                                                  Real & Equistrain_rate,
                                                                                  Real & Kapaq,
                                                                                  Real dt,
                                                                                  Real & crit_strain,
                                                                                  Real & crit_strain_rate,
                                                                                  Real & rdr_damage,
                                                                                  Real & nb_damage) {
  // __attribute__((unused)) Real dt
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
  // Real cst11 = 0.04;
   //Real cst11 = Alpha*((1.-0.5)*exp(-pow(Equistrain_rate*(1/8e5),10))+0.5);
  // Real cst11 = Alpha*(1.+(Equistrain_rate/Beta)); 
  //  Real cst11 = 0;
  //  if (Equistrain_rate >= 750000.){
  // cst11 = Alpha*(1.+(Equistrain_rate/500000.)); 
  //} else {
  // cst11 = Alpha*(1.+(Equistrain_rate/1000000.)); 
  //}
 
  // epsilono
  //     Real cst10 = 75./2000.;
         Real cst10 = Kapa0;
  //     Real cst10 = (Kapa0)*(1.+(Equistrain_rate/Beta)); 


  // Real cst10 = (75./3090.)*(1.+(std::abs(Equistrain_rate)/Beta))*(1+pow(std::abs(Equistrain)/0.1,2)); 
  //Real cst10 = (75./3090.)*(1.+(std::abs(Equistrain_rate)/Beta))*(0.17/std::max(std::abs(Equistrain),0.085)); 
  //  Real cst10 = (75./3090.)*(1.+(std::abs(Equistrain_rate)/Beta)); 
  // Real cst10 = (75./3090.)*(1.+(std::abs(Equistrain_rate)/Beta))*((1-dam)*0.17/std::max(std::abs(Equistrain),0.085));
  // Real cst10 = (75./3090.)*(1.+(std::abs(Equistrain_rate)/Beta))*((1-dam)*0.17/((0.085+(1+tanh(10000*(std::abs(Equistrain)-0.085)))/2)*(std::abs(Equistrain)-0.085)));

  //    Real cst10 = Equistrain_rate; 
  //Real cst10 = (75./3090. - Alpha + 0.01)*exp(-pow(Equistrain_rate/(5e-5),10))+ Alpha - 0.01;
  //test
  //       Real sigmaco=(this->E)*Kapa0*Alpha/(Alpha-Kapa0);
  //       Real sigmaco=(this->E)*(75./3090.)*Alpha/(Alpha-75./3090);

  //    Real cst10= (sigmaco/this->E)*(1.+(Equistrain_rate/Beta))/(1.+sigmaco/(Alpha*this->E));
  //deg2
  //    Real cst10 = (75./3090.)*(1.+(Equistrain_rate/Beta)+(Equistrain_rate/(Beta*2))*(Equistrain_rate/(Beta*2))); 
  //ln

  //new law

    


  // Real cst10 = (Kapa0)*(1.+(log(Equistrain_rate+1.)/(Beta/10000.))); 

  // epsilonc
       Real cst11 = Alpha;
       //Real cst11 = 1000*Alpha*(1.-(Equistrain_rate/Beta)); 

//---------------------------------------------------------------------------------------
//
//   Real cst9 = 0;
   // if (cst10 > Kapa0) {
  if (cst10 > Kapa0) {
    // if (cst10 < 10.) {
      if (cst10 < 0.16) {
	//cst9=cst10;
    } else {
	//cst9=10.;
	//cst9=0.16;
    }
  } else {
    //    cst9 = Kapa0;
  }
//
//  Real Fd1 = Equistrain - Kapa0;
  Real Fd1 = Equistrain - Kapaq;
    if (Fd1 > 0)
     {
//       //Real Fd = Equistrain - Kapaq;
//       // if (Fd  > 0)
//       //{
//	      Real dam1 = 1. - cst10/Equistrain * ((cst11-Equistrain)/(cst11 - cst10));
	     Real dam1 = 1. - Kapaq/Equistrain * ((cst11-Equistrain)/(cst11 - Kapaq));
	     //	     Real dam1 = 1. - cst9/Equistrain * ((cst11-Equistrain)/(cst11 - cst9));
        if (dam1 > dam){
	  dam = std::min(dam1,1.);
           if (rdr_damage < 1000.){
             nb_damage = nb_damage + 1.;
	   }
	   if (dam >= 1. && rdr_damage < 1000.){
	     	crit_strain = Equistrain;
	     	crit_strain_rate = Equistrain_rate;
              rdr_damage = 2000.;
	   }
//	     //           if (dam1 >= 0.9999){
//  	     //dam =1.;
//	     //	   }
//	     //else {    
//  	     //dam = std::min(dam1,1.);
//             //crit_strain=cst9;
//	     //}
//	     //	     }
	}
     }
//---------------------------------------------------------------------------------------
// last model
//
//---------------------------------------------------------------------------------------
//
//    Real epsilonoi = 75./3090.;
//    Real epsilonc = 0.17;
//    Real epsilonmax = 0.10;
//    Real deltaepsilonc = 0.03;
//    Real sE = 3090e6;
//    
//    Real epsilononew = epsilonoi * ( 1 + (Equistrain_rate / Beta));
//    Real epsilono = 0.;
//
//    if (epsilononew > epsilonoi) {
//      if (epsilononew < epsilonmax - 0.01){
//	epsilono = epsilononew;
//      } else {
//        epsilono = epsilonmax - 0.01;
//      }
//    } else {
//      epsilono = epsilonoi;
//    }
//
//    Real fepsilonmax = sE*epsilonmax / ( epsilonmax - epsilonc );
//    Real fepsilonoi = sE*epsilonoi / ( epsilonoi - epsilonc );
//
//    Real epsilonm = ( fepsilonmax*epsilonc + epsilono*( sE - fepsilonoi ) ) / ( fepsilonmax - fepsilonoi );
//    Real epsilona = ( epsilonmax * ( epsilonm - epsilonc ) ) / ( epsilonmax - epsilonc );
//
//    //----------------
//
//    Real Fd1 = Equistrain - epsilono; 
//    if (Fd1 > 0)
//     {
//       Real Fd2 = Equistrain - epsilonm;
//       if (Fd2 > 0)
//	 {
//	   Real dam3 = 1 - (epsilona / Equistrain) + (epsilonmax * (Equistrain - epsilonm)) / (deltaepsilonc * Equistrain); 
//           if (dam3 > dam){
//   	     dam = std::min(dam3,1.);
//             crit_strain=epsilonm;
//  	   }
//	 } else {
//	 Real dam2 = 1 - ((epsilona - epsilono) * (Equistrain - epsilono)) / (( epsilonm - epsilono) * Equistrain) - epsilono / Equistrain;
//           if (dam2 > dam){
//   	     dam = std::min(dam2,1.);
//             crit_strain=epsilonm;
//	   }
//       }
//     }
//
//---------------------------------------------------------------------------------------
// new model realistic
//---------------------------------------------------------------------------------------
//
//  // Real sgamma = 0.04;
//  Real sn1 = 1.;
//  Real sn2 = 0.5;
//  Real betal = 600000.;
//   // Real sn = 10.;
//  // Real betal = 5000000.;
//  Real epsilono = 75./3090.;
//  Real sE = 3090e6;
//  //    Real sbeta = 0.02;
//  // ok this is work
//  // Real epsilonmax = 0.05;
//  // Real epsilonco = 0.055;
//
//  // it is the best
//  //Real epsilonmax = 0.0525;
//  //Real epsilonco = 0.0575;
//
//  Real epsilonmax = 0.0525;
//  Real epsilonco = 0.0575;
//
//
//
//
//  Real sigmaco = 78e6;
//  Real deltaepsilonc = 0.0001;
//  // Real deltaepsilonc = 0.005178;
//
//  //   Real epsilonmnew = sbeta*pow(1/Equistrain_rate,sn)+sgamma;
//
//   Real epsilonmnew = (epsilonco-epsilonmax)*exp(-pow(Equistrain_rate*(1/betal),sn2))+epsilonmax;
//   Real epsilononew = (epsilono-epsilonmax)*exp(-pow(Equistrain_rate*(1/betal),sn1))+epsilonmax;
//  // Real epsilonmnew = epsilonco;
//  // Real epsilononew = epsilono;
//
//  Real epsilonm = 0.;
//  Real epsilonoc = 0.;
//
//  if (epsilonmnew > epsilonmax) {
//    if (epsilonmnew < epsilonco){
//	epsilonm = epsilonmnew;
//    } else {
//      epsilonm = epsilonco;
//    }
//  } else {
//    epsilonm = epsilonmax;
//  }
//
//  if (epsilononew > epsilono) {
//    if (epsilononew < epsilonmax){
//	epsilonoc = epsilonmnew;
//    } else {
//      epsilonoc = epsilonmax;
//    }
//  } else {
//    epsilonoc = epsilono;
//  }
//
//
//  Real epsiloncoprim = sigmaco/sE;
//  Real epsilona = epsilonmax + (epsiloncoprim-epsilonmax)*(epsilonm-epsilonmax)/(epsilonco-epsilonmax);
//  Real sigmam= sE*epsilona;
//  Real epsilonc = sigmam*deltaepsilonc/sigmaco+epsilonm;
//
//  //----------------
//
//  Real Fd1 = Equistrain - epsilono; 
//  if (Fd1 > 0)
//   {
//     Real Fd2 = Equistrain - epsilonm;
//     if (Fd2 > 0)
//	 {
//	   Real dam3 = 1-(epsilona/Equistrain)*(1-(Equistrain-epsilonm)/(epsilonc-epsilonm));
//         if (dam3 > dam){
// 	     dam = std::min(dam3,1.);
//	      if (dam < 0.9999){
//		crit_strain = std::max(Equistrain_rate,crit_strain);
//	      }
//	     }
//	 } else {
//	 Real dam2 = 1-(epsilonoc/Equistrain)-(Equistrain-epsilonoc)/(epsilonm-epsilonoc)*(epsilona-epsilonoc)/Equistrain;
//         if (dam2 > dam){
// 	     dam = std::min(dam2,1.);
//	      if (dam < 0.9999){
//		crit_strain = std::max(Equistrain_rate,crit_strain);
//	      }
//	    }
//	   }
//   }
//
//
//---------------------------------------------------------------------------------------
// publication Wu
//---------------------------------------------------------------------------------------
//  
//  // Real sgamma = 0.04;
////Real sn1 = 0.2;
//	       Real sn2 = 0.1;
//    //Real betal1 = 1e2;
//    Real betal2 = 100.;
//	//  Real sn1 = 2.;
//	//        Real sn1 = 0.5;
//	// Real sn2 = 2.;
//	// Real betal1 = 6e5;
//	//Real betal2 = 6e5;
//
//   // Real sn = 10.;
//  // Real betal = 5000000.;
//  Real epsilono = 1./4000.;
//  Real sE = 4000e6;
//  //    Real sbeta = 0.02;
//  // ok this is work
//  // Real epsilonmax = 0.04;
//  // Real epsilonco = 0.1;
//
//
//  Real epsilonmax = 0.03;
//  Real epsilonco = 0.055;
//
//
//  Real sigmaco = 75e6;
//  Real deltaepsilonc = 0.04;
//  // Real deltaepsilonc = 0.002045;
//
//  //   Real epsilonmnew = sbeta*pow(1/Equistrain_rate,sn)+sgamma;
//
//   Real epsilonmnew = (epsilonco-epsilonmax)*exp(-pow(Equistrain_rate*(1/betal2),sn2))+epsilonmax;
//   //  Real epsilononew = (epsilonmax-epsilonco)*exp(-pow(Equistrain_rate*(1/betal1),sn1))+epsilonco;
//  // Real epsilonmnew = epsilonco;
//   Real epsilononew = epsilono;
//
//  Real epsilonm = 0.;
//  Real epsilonoc = 0.;
//
//  if (epsilonmnew > epsilonmax) {
//    if (epsilonmnew < epsilonco){
//	epsilonm = epsilonmnew;
//    } else {
//      epsilonm = epsilonco;
//    }
//  } else {
//    epsilonm = epsilonmax;
//  }
//
//  if (epsilononew > epsilono) {
//    if (epsilononew < epsilonmax){
//	epsilonoc = epsilonmnew;
//    } else {
//      epsilonoc = epsilonmax;
//    }
//  } else {
//    epsilonoc = epsilono;
//  }
//
//
//  Real epsiloncoprim = sigmaco/sE;
//  Real epsilona = epsilonmax + (epsiloncoprim-epsilonmax)*(epsilonm-epsilonmax)/(epsilonco-epsilonmax);
//  Real sigmam= sE*epsilona;
//  Real epsilonc = sigmam*deltaepsilonc/sigmaco+epsilonm;
//
//  //----------------
//
//  Real Fd1 = Equistrain - epsilono; 
//  if (Fd1 > 0)
//   {
//     Real Fd2 = Equistrain - epsilonm;
//     if (Fd2 > 0)
//	 {
//	   Real dam3 = 1-(epsilona/Equistrain)*(1-(Equistrain-epsilonm)/(epsilonc-epsilonm));
//         if (dam3 > dam){
// 	     dam = std::min(dam3,1.);
//       if (rdr_damage < 1000.){
//         nb_damage = nb_damage + 1.;
//	       }
//	     // 	      if (dam < 0.9999){
//	     //	crit_strain = std::max(Equistrain,crit_strain);
//	     //	crit_strain_rate = std::max(Equistrain_rate,crit_strain_rate);
//	     if (dam >= 1. && rdr_damage<1000.){
//	     	crit_strain = Equistrain;
//	     	crit_strain_rate = Equistrain_rate;
//          rdr_damage = 2000.;
//		}
//	      }
//	 } else {
//	 Real dam2 = 1-(epsilonoc/Equistrain)-(Equistrain-epsilonoc)/(epsilonm-epsilonoc)*(epsilona-epsilonoc)/Equistrain;
//         if (dam2 > dam){
// 	     dam = std::min(dam2,1.);
//       if (rdr_damage < 1000.){
//         nb_damage = nb_damage + 1.;
//	       }
//	     // if (dam < 0.9999){
//	     //crit_strain = std::max(Equistrain,crit_strain);
//	     //	crit_strain_rate = std::max(Equistrain_rate,crit_strain_rate);
//	     if (dam >= 1. && rdr_damage < 1000.){
//	     	crit_strain = Equistrain;
//	     	crit_strain_rate = Equistrain_rate;
//          rdr_damage = 2000.;
//		 }
//	     }
//	    }
//	   }
//	 
//
//---------------------------------------------------------------------------------------
// good
//---------------------------------------------------------------------------------------
 
//---------------------------------------------------------------------------------------
// delay
//---------------------------------------------------------------------------------------
//    Real viscosity = 0.0001;
//    Real damrateinfini = 10000000000.;
//    Real gequi = 1. - cst10/Equistrain * ((cst11-Equistrain)/(cst11 - cst10));
//    Real damrate = damrateinfini * (1. - exp(-1*viscosity * (gequi - dam)));
//      if (damrate > 0) dam = dam + damrate*dt;
//
//---------------------------------------------------------------------------------------

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
