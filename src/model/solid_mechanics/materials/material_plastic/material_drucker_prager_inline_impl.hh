/**
 * @file   material_linear_isotropic_hardening_inline_impl.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Sep 09 2020
 * @date last modification: Wed Sep 09 2020
 *
 * @brief  Implementation of the inline functions of the material
 * Drucker-Prager
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "material_drucker_prager.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
/// Deviatoric Stress
template<UInt dim>  
inline void MaterialDruckerPrager<dim>::computeDeviatoricStress(const Matrix<Real> & sigma,
								Matrix<Real> & sigma_dev){
  for (UInt i = 0; i < dim; ++i){
    for (UInt j = 0; j < dim; ++j){
      sigma_dev(i, j) = sigma(i, j);
    }
  }

  sigma_dev -= Matrix<Real>::eye(dim, sigma.trace() / dim);
}

/* -------------------------------------------------------------------------- */
/// Yield Stress
template<UInt dim>
inline Real MaterialDruckerPrager<dim>::computeYieldStress(const Matrix<Real> &
							   sigma) {
  return this->alpha * sigma.trace() - this->k;
}


/* -------------------------------------------------------------------------- */  
/// Yield function
template<UInt dim>
inline Real MaterialDruckerPrager<dim>::computeYieldFunction(const Matrix<Real> &
							     sigma) {
  Matrix<Real> sigma_dev(dim, dim, 0);
  this->computeDeviatoricStress(sigma, sigma_dev);
  
  // compute deviatoric invariant J2
  Real j2 = (1./ 2.) * sigma_dev.doubleDot(sigma_dev);
  Real sigma_dev_eff = std::sqrt(3. * j2);

  Real modified_yield_stress = computeYieldStress(sigma);

  return sigma_dev_eff + modified_yield_stress;
}

/* -------------------------------------------------------------------------- */
/// Threshold Stress  
template<UInt dim>
inline bool MaterialDruckerPrager<dim>::aboveThresholdStress(const Matrix<Real> & sigma,
							     Real & k, Real & alpha) {
  Matrix<Real> sigma_dev(dim, dim, 0);
  this->computeDeviatoricStress(sigma, sigma_dev);

  Real I1 = sigma.trace();
  Real J2 = (1./2) *sigma_dev.doubleDot(sigma_dev);
    
  if(I1 > k/alpha) {
    Real beta = sqrt(3.*J2)/(I1-k/alpha);
    Real gamma = tan(M_PI/2. - atan(alpha));
      
    if(beta < gamma)
      return true;
  }
    
  return false;    
}

/* -------------------------------------------------------------------------- */
template<UInt dim>  
inline Real MaterialDruckerPrager<dim>::computeObjectiveFunction(const Matrix<Real> & sigma_guess,
								 const Matrix<Real> & sigma_trial,
								 Real & plastic_multiplier_guess,
								 Vector<Real> & gradient_f,
								 Vector<Real> & delta_inelastic_strain,
								 Real & yield_function){

  UInt size = voigt_h::size;
  
  // elastic stifnness tensor
  Matrix<Real> De(size, size, 0.);
  MaterialElastic<dim>::computeTangentModuliOnQuad(De);

  // elastic compliance tensor
  Matrix<Real> Ce(size, size, 0.);
  Ce.inverse(De);

  // objective function to be computed
  Vector<Real> f(size, 0.);

  // compute gradient
  computeJacobian(sigma_guess, gradient_f);
  
  // compute yield function  
  yield_function =  computeYieldFunction(sigma_guess);
    
  // compute increment strain
  auto sigma_trial_voigt = voigt_h::matrixToVoigt(sigma_trial);
  auto sigma_guess_voigt = voigt_h::matrixToVoigt(sigma_guess);
  auto tmp = sigma_trial_voigt - sigma_guess_voigt;
  delta_inelastic_strain.mul<false>(Ce, tmp);

  // compute objective function 
  f.mul<false>(De, gradient_f, plastic_multiplier_guess);
  f = tmp - f; 
    
  // compute error
  Real error = std::max(f.norm<L_2>(), std::abs(yield_function));
  return error;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialDruckerPrager<dim>::computeJacobian(const Matrix<Real> & sigma_guess,
							Vector<Real> & gradient_f) {

  UInt size = voigt_h::size;

  // scaling matrix for computing gradient and hessian from voigt notation
  Matrix<Real> scaling_matrix(size, size, 0.);
  scaling_matrix.eye(1.);
  for(auto i : arange(dim, size)) {
    for(auto j : arange(dim, size)){
      scaling_matrix(i, j) *= 2.;
    }
  }

  // krocker delta vector in voigt notation
  Vector<Real> kronecker_delta(size, 0.);
  for(auto i : arange(dim)) {
    kronecker_delta[i] = 1.;
  }

  // compute deviatoric invariant
  Matrix<Real> sigma_dev(dim, dim, 0);
  
  this->computeDeviatoricStress(sigma_guess, sigma_dev);
  
  Real j2 = (1./2.) * sigma_dev.doubleDot(sigma_dev);

  Vector<Real> sigma_dev_voigt = voigt_h::matrixToVoigt(sigma_dev);

  gradient_f.mul<false>(scaling_matrix, sigma_dev_voigt, 3./ (2. * std::sqrt(3. * j2)) );
  gradient_f += alpha * kronecker_delta;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>  
inline void MaterialDruckerPrager<dim>::computeHessian(const Matrix<Real> & sigma_guess,
						       Matrix<Real> & hessian_f) {

  UInt size = voigt_h::size;

  // scaling matrix for computing gradient and hessian from voigt notation
  Matrix<Real> scaling_matrix(size, size, 0.);
  scaling_matrix.eye(1.);
  for(auto i : arange(dim, size)) {
    for(auto j : arange(dim, size)){
      scaling_matrix(i, j) *= 2.;
    }
  }

  // krocker delta vector in voigt notation
  Vector<Real> kronecker_delta(size, 0.);
  for(auto i : arange(dim)) {
    kronecker_delta[i] = 1.;
  }

  // compute deviatoric invariant J2
  Matrix<Real> sigma_dev(dim, dim, 0);

  this->computeDeviatoricStress(sigma_guess, sigma_dev);
   
  Real j2 = (1./2.) * sigma_dev.doubleDot(sigma_dev);
  
  auto sigma_dev_voigt = voigt_h::matrixToVoigt(sigma_dev);

  Vector<Real> temp(sigma_dev_voigt.size());
  temp.mul<false>(scaling_matrix, sigma_dev_voigt);

  Matrix<Real> id(kronecker_delta.size(), kronecker_delta.size());
  id.outerProduct(kronecker_delta, kronecker_delta);
  id *=  -1./3.;
  id += Matrix<Real>::eye(kronecker_delta.size(), 1.);

  Matrix<Real> tmp3(kronecker_delta.size(), kronecker_delta.size());
  tmp3.mul<false, false>(scaling_matrix, id);
  hessian_f.outerProduct(temp, temp);
  hessian_f *= -9./(4.* pow(3.*j2, 3./2.));
  hessian_f += (3./(2.* pow(3.*j2, 1./2.)))*tmp3;
}

/* -------------------------------------------------------------------------- */ 
template<UInt dim>
inline void MaterialDruckerPrager<dim>::computeGradientAndPlasticMultplier(
       const Matrix<Real> & sigma_trial, Real & plastic_multiplier_guess,      
       Vector<Real> & gradient_f, Vector<Real> & delta_inelastic_strain,
       UInt max_iterations, 
       Real tolerance) {

  UInt size = voigt_h::size;

  // guess stress state at each iteration, initial guess is the trial state
  Matrix<Real> sigma_guess(sigma_trial);

  // hessian matrix of yield surface
  Matrix<Real> hessian_f(size, size, 0.);

  // variation in sigma at each iteration
  Vector<Real> delta_sigma(size, 0.);
  
  // plastic multiplier guess at each iteration, initial guess is zero 
  plastic_multiplier_guess = 0.;
  
  // gradient of yield surface in voigt notation
  gradient_f.zero();
  
  // plastic strain increment at each iteration 
  delta_inelastic_strain.zero();

  // elastic stifnness tensor
  Matrix<Real> De(size, size, 0.);
  MaterialElastic<dim>::computeTangentModuliOnQuad(De);

  // elastic compliance tensor
  Matrix<Real> Ce(size, size, 0.);
  Ce.inverse(De);
     
  // yield function value at each iteration
  Real yield_function;

  /* ------------------------------------------------------- */
  /* Projecting stress state to origin of yield function     */
  /* if first invariant is greater than the threshold        */
  /* ------------------------------------------------------- */

  if(this->aboveThresholdStress(sigma_guess, k, alpha) and
     this->alpha > 0) {

    this->projectionOnThreshold(sigma_guess, gradient_f,
				plastic_multiplier_guess, tolerance);
   
    auto delta_sigma_final = sigma_trial - sigma_guess;
    auto delta_sigma_voigt = voigt_h::matrixToVoigt(delta_sigma_final);
    
    delta_inelastic_strain.mul<false>(Ce, delta_sigma_voigt);
    return;
  }
    
  /* -------------------------------------------------- */
  /* If Drucker-prager is Capped                        */
  /* -------------------------------------------------- */

  /*Real I1 = sigma_guess.trace();
  
  if(is_capped) {

    Matrix<Real> sigma_dev(dim, dim, 0);
    this->computeDeviatoricStress(sigma_guess, sigma_dev);
    Real J2 = (1./2) *sigma_dev.doubleDot(sigma_dev);

    
    if(I1 < hydrostatic_min and sqrt(3.*J2) < deviatoric_max) {

      auto projection_error = I1 - hydrostatic_min;

      while(tolerance < std::abs(projection_error)) {

	Matrix<Real> delta_sigma(dim, dim);
	Matrix<Real> jacobian(dim, dim);
	Matrix<Real> jacobian_inv(dim, dim);

	jacobian = Matrix<Real>::eye(dim, dim);
	jacobian_inv.inverse(jacobian);

	delta_sigma += -projection_error * jacobian_inv;
	sigma_guess += delta_sigma;
	projection_error = sigma_guess.trace() - hydrostatic_min;
      }

      auto delta_sigma_final = sigma_trial - sigma_guess;
      auto delta_sigma_voigt = voigt_h::matrixToVoigt(delta_sigma_final);
    
      delta_inelastic_strain.mul<false>(Ce, delta_sigma_voigt);
      return;
    }
        

    Real beta = (sqrt(3.*J2) - deviatoric_max )/(I1- hydrostatic_min);
    Real gamma = tan(M_PI/2. - atan(alpha));

    
    if((I1 > hydrostatic_min and beta > gamma) or
       (I1 < hydrostatic_min and sqrt(3.*J2) > deviatoric_max )) {

      // lambda function to compute gradient of yield surface in voigt notation
      auto compute_gradient_f3 = [&sigma_guess, &scaling_matrix, &kronecker_delta,
				  &gradient_f](){
	
	Matrix<Real> sigma_dev(dim, dim, 0);

	for (UInt i = 0; i < dim; ++i)
	  for (UInt j = 0; j < dim; ++j)
	    sigma_dev(i, j) = sigma_guess(i, j);

	sigma_dev -= Matrix<Real>::eye(dim, sigma_guess.trace() / dim);
    
	Vector<Real> sigma_dev_voigt = voigt_h::matrixToVoigt(sigma_dev);

	// compute deviatoric invariant
	Real j2 = (1./2.) * sigma_dev.doubleDot(sigma_dev);
        
	gradient_f.mul<false>(scaling_matrix, sigma_dev_voigt, 3./ (2. * std::sqrt(3. * j2)) );
      
      };

      auto projection_error = sqrt(3.*J2) - deviatoric_max;

      while( tolerance < projection_error ) {
	compute_gradient_f3();

	Vector<Real> tmp1(size);
	tmp1.mul<false>(Ce, gradient_f);

	auto denominator = gradient_f.dot(tmp1);

	auto delta_plastic_multiplier = projection_error;
	delta_plastic_multiplier /= denominator;
    
	plastic_multiplier_guess += delta_plastic_multiplier;
	delta_sigma.mul<false>(Ce, gradient_f, delta_plastic_multiplier); 
    
	// compute sigma_guess from voigt notation
	Matrix<Real> delta_sigma_mat(dim, dim);
	voigt_h::voigtToMatrix(delta_sigma, delta_sigma_mat);
	sigma_guess -= delta_sigma_mat;

	this->computeDeviatoricStress(sigma_guess, sigma_dev);
	J2 = (1./2) *sigma_dev.doubleDot(sigma_dev);

	projection_error = sqrt(3.*J2) - deviatoric_max;
      }

      projection_error = sigma_guess.trace() - hydrostatic_min;

      while(tolerance < std::abs(projection_error)) {

	Matrix<Real> delta_sigma(dim, dim);
	Matrix<Real> jacobian(dim, dim);
	Matrix<Real> jacobian_inv(dim, dim);

	jacobian = Matrix<Real>::eye(dim, dim);
	jacobian_inv.inverse(jacobian);

	delta_sigma += -projection_error * jacobian_inv;
	sigma_guess += delta_sigma;
	projection_error = sigma_guess.trace() - hydrostatic_min;
      }


      auto delta_sigma_final = sigma_trial - sigma_guess;
      auto delta_sigma_voigt = voigt_h::matrixToVoigt(delta_sigma_final);
    
      delta_inelastic_strain.mul<false>(Ce, delta_sigma_voigt);

      return;
    }
    
  }*/
  
  /* -------------------------------------------------- */
  /* Generalized cutting plane algorithm (explicit)     */
  /* -------------------------------------------------- */
  /*UInt iteration = 0;
  
  auto error = this->computeYieldFunction(sigma_guess);

  while( tolerance < error ) {
    compute_gradient_f(alpha);

    Vector<Real> tmp1(size);
    tmp1.mul<false>(Ce, gradient_f);

    auto denominator = gradient_f.dot(tmp1);

    auto delta_plastic_multiplier = error;
    delta_plastic_multiplier /= denominator;
    
    plastic_multiplier_guess += delta_plastic_multiplier;
    delta_sigma.mul<false>(Ce, gradient_f, delta_plastic_multiplier); 
    
    // compute sigma_guess from voigt notation
    Matrix<Real> delta_sigma_mat(dim, dim);
    voigt_h::voigtToMatrix(delta_sigma, delta_sigma_mat);
    sigma_guess -= delta_sigma_mat;

    // compute increment inelastic strain
    auto sigma_trial_voigt = voigt_h::matrixToVoigt(sigma_trial);
    auto sigma_guess_voigt = voigt_h::matrixToVoigt(sigma_guess);
    auto tmp = sigma_trial_voigt - sigma_guess_voigt;
    delta_inelastic_strain.mul<false>(Ce, tmp);

    
    error = this->computeYieldFunction(sigma_guess);
    
    iteration++;
  }*/

  /* ----------------------------------------------------------------------- */
  /* Projection on yield surface using Closet point projection (implicit)    */
  /* ----------------------------------------------------------------------- */ 

  Real projection_error = computeObjectiveFunction(sigma_guess, sigma_trial, plastic_multiplier_guess,
  						   gradient_f, delta_inelastic_strain, yield_function);
  
  /* --------------------------- */
  /* iteration loop              */
  /* --------------------------- */

  Matrix<Real> xi(size, size);

  Matrix<Real> xi_inv(size, size);

  Vector<Real> tmp(size);

  Vector<Real> tmp1(size);

  Matrix<Real> tmp2(size, size);
  
  UInt iterations{0};
  while(tolerance < projection_error and iterations < max_iterations) {
    
    // compute hessian at previous step
    computeHessian(sigma_guess, hessian_f);
    
    // compute inverse matrix Xi    
    xi = Ce + plastic_multiplier_guess * hessian_f;

    // compute inverse matrix Xi
    xi_inv.inverse(xi);
    
    tmp.mul<false>(xi_inv, gradient_f);

    auto denominator = gradient_f.dot(tmp);
    
    // compute plastic multplier guess
    tmp1.mul<false>(xi_inv, delta_inelastic_strain);
    plastic_multiplier_guess = gradient_f.dot(tmp1);
    plastic_multiplier_guess += yield_function;
    plastic_multiplier_guess /= denominator;

    // compute delta sigma    
    tmp2.outerProduct(tmp, tmp);
    tmp2 /= denominator;

    tmp2 = xi_inv - tmp2;
    delta_sigma.mul<false>(tmp2, delta_inelastic_strain);

    delta_sigma -= tmp*yield_function/denominator;

    // compute sigma_guess
    Matrix<Real> delta_sigma_mat(dim, dim);
    voigt_h::voigtToMatrix(delta_sigma, delta_sigma_mat);
    sigma_guess += delta_sigma_mat;

    projection_error = computeObjectiveFunction(sigma_guess, sigma_trial, plastic_multiplier_guess,
						gradient_f, delta_inelastic_strain, yield_function);

    iterations++;
  }

  if ( iterations >= max_iterations and tolerance < projection_error) {
    AKANTU_ERROR("Maximum iterations reached in projection without reaching the tolerance "
		 << " Increase the iterations or reduce the tolerance "
		 << " Error is " << projection_error);
  }

}

/* -------------------------------------------------------------------------- */  
/// Infinitesimal deformations
template <UInt dim>
inline void MaterialDruckerPrager<dim>::computeStressOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
    Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
    Matrix<Real> & inelastic_strain,  const Matrix<Real> & previous_inelastic_strain,
    const Real & sigma_th, const Real & previous_sigma_th) {
  
  Real delta_sigma_th = sigma_th - previous_sigma_th;

  Matrix<Real> grad_delta_u(grad_u);
  grad_delta_u -= previous_grad_u;

  // Compute trial stress, sigma_tr
  Matrix<Real> sigma_tr(dim, dim);
  MaterialElastic<dim>::computeStressOnQuad(grad_delta_u, sigma_tr,
                                            delta_sigma_th);
  sigma_tr += previous_sigma;

  // Compute the yielding sress  
  bool initial_yielding =
    (this->computeYieldFunction(sigma_tr) > 0);

  bool above_capped =
    is_capped and (sigma_tr.trace() < hydrostatic_min);

  // use closet point projection to compute the plastic multiplier and
  // gradient and inealstic strain at the surface for the given trial
  // stress state
  Matrix<Real> delta_inelastic_strain(dim, dim, 0.);

  if(initial_yielding or above_capped) {

    UInt size = voigt_h::size;

    // plastic multiplier
    Real dp{0.};

    // gradient of yield surface in voigt notation
    Vector<Real> gradient_f(size, 0.);

    // inelastic strain in voigt notation
    Vector<Real> delta_inelastic_strain_voigt(size, 0.);

    // compute using closet-point projection
    this->computeGradientAndPlasticMultplier(sigma_tr, dp, gradient_f,
					     delta_inelastic_strain_voigt);

    for(auto i: arange(dim, size))
      delta_inelastic_strain_voigt[i] /= 2.;
    
    voigt_h::voigtToMatrix(delta_inelastic_strain_voigt,
			   delta_inelastic_strain);
  }
  

  // Compute the increment in inelastic strain
  MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(
      grad_delta_u, sigma, previous_sigma, inelastic_strain,
      previous_inelastic_strain, delta_inelastic_strain);
}

/* -------------------------------------------------------------------------- */  
template<UInt dim>
inline void MaterialDruckerPrager<dim>::projectionOnThreshold(Matrix<Real> & sigma_guess,
								 Vector<Real> & gradient_f,
								 Real & plastic_multiplier_guess,
								 Real & tolerance){

  UInt size = voigt_h::size;

  // elastic stifnness tensor
  Matrix<Real> De(size, size, 0.);
  MaterialElastic<dim>::computeTangentModuliOnQuad(De);

  // elastic compliance tensor
  Matrix<Real> Ce(size, size, 0.);
  Ce.inverse(De);

  
  // variation in sigma at each iteration
  Vector<Real> delta_sigma(size, 0.);

  
  // krocker delta vector in voigt notation
  Vector<Real> kronecker_delta(size, 0.);
  for(auto i : arange(dim)) {
    kronecker_delta[i] = 1.;
  }
  
  // hessian matrix of yield surface
  Matrix<Real> hessian_f(size, size, 0.);

  // scaling matrix for computing gradient and hessian from voigt notation
  Matrix<Real> scaling_matrix(size, size, 0.);
  scaling_matrix.eye(1.);
  for(auto i : arange(dim, size)) {
    for(auto j : arange(dim, size)){
      scaling_matrix(i, j) *= 2.;
    }
  }
 

  
  auto update_first_obj = [&sigma_guess]() {
      
      Matrix<Real> sigma_dev(dim, dim, 0);

      for (UInt i = 0; i < dim; ++i)
	for (UInt j = 0; j < dim; ++j)
	  sigma_dev(i, j) = sigma_guess(i, j);

      sigma_dev -= Matrix<Real>::eye(dim, sigma_guess.trace() / dim);

      Real J2 = (1./2) *sigma_dev.doubleDot(sigma_dev);
      
      auto error = sqrt(3.*J2);
      
      return error;
    };

    auto update_sec_obj = [&sigma_guess](Real & k, Real & alpha) {
      auto error = alpha*sigma_guess.trace() - k;
      return error;
    };


    // lambda function to compute gradient of yield surface in voigt notation
    auto compute_gradient_f2 = [&sigma_guess, &scaling_matrix, &kronecker_delta,
				&gradient_f](){
    
      Matrix<Real> sigma_dev(dim, dim, 0);

      for (UInt i = 0; i < dim; ++i)
	for (UInt j = 0; j < dim; ++j)
	  sigma_dev(i, j) = sigma_guess(i, j);

      sigma_dev -= Matrix<Real>::eye(dim, sigma_guess.trace() / dim);
    
      Vector<Real> sigma_dev_voigt = voigt_h::matrixToVoigt(sigma_dev);

      // compute deviatoric invariant
      Real j2 = (1./2.) * sigma_dev.doubleDot(sigma_dev);
        
      gradient_f.mul<false>(scaling_matrix, sigma_dev_voigt, 3./ (2. * std::sqrt(3. * j2)) );
      
    };

    auto projection_error = update_first_obj();

    while( tolerance < projection_error ) {
      compute_gradient_f2();

      Vector<Real> tmp1(size);
      tmp1.mul<false>(Ce, gradient_f);

      auto denominator = gradient_f.dot(tmp1);

      auto delta_plastic_multiplier = projection_error;
      delta_plastic_multiplier /= denominator;
    
      plastic_multiplier_guess += delta_plastic_multiplier;
      delta_sigma.mul<false>(Ce, gradient_f, delta_plastic_multiplier); 
    
      // compute sigma_guess from voigt notation
      Matrix<Real> delta_sigma_mat(dim, dim);
      voigt_h::voigtToMatrix(delta_sigma, delta_sigma_mat);
      sigma_guess -= delta_sigma_mat;

      projection_error = update_first_obj();
    }


    projection_error = update_sec_obj(k, alpha);

    while(tolerance < projection_error) {

      Matrix<Real> delta_sigma(dim, dim);
      Matrix<Real> jacobian(dim, dim);
      Matrix<Real> jacobian_inv(dim, dim);

      jacobian = this->alpha * Matrix<Real>::eye(dim, dim);
      jacobian_inv.inverse(jacobian);

      delta_sigma += -projection_error * jacobian_inv;
      sigma_guess += delta_sigma;
      projection_error = update_sec_obj(k, alpha);
    } 
}
/* -------------------------------------------------------------------------- */
/// Finite deformations
template <UInt dim>
inline void MaterialDruckerPrager<dim>::computeStressOnQuad(
    __attribute__((unused)) const Matrix<Real> & grad_u,
    __attribute__((unused)) const Matrix<Real> & previous_grad_u,
    __attribute__((unused)) Matrix<Real> & sigma,
    __attribute__((unused)) const Matrix<Real> & previous_sigma,
    __attribute__((unused)) Matrix<Real> & inelastic_strain,
    __attribute__((unused)) const Matrix<Real> & previous_inelastic_strain,
    __attribute__((unused)) const Real & sigma_th,
    __attribute__((unused)) const Real & previous_sigma_th,
    __attribute__((unused)) const Matrix<Real> & F_tensor) {

  
}

/* -------------------------------------------------------------------------- */

template <UInt dim>
inline void MaterialDruckerPrager<dim>::computeTangentModuliOnQuad(
    __attribute__((unused)) Matrix<Real> & tangent,
    __attribute__((unused)) const Matrix<Real> & grad_u,
    __attribute__((unused)) const Matrix<Real> & previous_grad_u,
    __attribute__((unused)) const Matrix<Real> & sigma_tensor,
    __attribute__((unused)) const Matrix<Real> & previous_sigma_tensor) const {

}  
  
}
