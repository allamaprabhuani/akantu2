/**
 * @file   material_cohesive_Exponential.cc
 * @author Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Feb 20 12:14:16 2012
 *
 * @brief  Exponential irreversible cohesive law of mixed mode loading
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
#include "material_cohesive_exponential.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveExponential<spatial_dimension>::MaterialCohesiveExponential(SolidMechanicsModel & model, const ID & id) :
  MaterialCohesive(model,id) {
  AKANTU_DEBUG_IN();

  sigma_c   = 0;
  beta      = 0;
  G_cI      = 0;
  G_cII     = 0;

  // this->initInternalVector(delta_max, 1, _ek_cohesive);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveExponential<spatial_dimension>::~MaterialCohesiveExponential() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  kappa   = G_cII / G_cI;
  delta_c = 2 * G_cI / sigma_c;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::resizeCohesiveVectors() {
  MaterialCohesive::resizeCohesiveVectors();
  //this->resizeInternalVector(delta_max, _ek_cohesive);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialCohesiveExponential<spatial_dimension>::setParam(const std::string & key, 
				      const std::string & value,
				      const ID & id) {
  std::stringstream sstr(value);
  if(key == "sigma_c") { sstr >> sigma_c; }
  else if(key == "beta") { sstr >> beta; }
  else if(key == "G_cI") { sstr >> G_cI; }
  else if(key == "G_cII") { sstr >> G_cII; }
  else { return Material::setParam(key, value, id); }
  return true;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::computeTraction(const Vector<Real> & normal,
						  ElementType el_type,
						  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  Vector<Real>::iterator<types::RVector> traction_it =
    tractions(el_type, ghost_type).begin(spatial_dimension);

  Vector<Real>::iterator<types::RVector> opening_it =
    opening(el_type, ghost_type).begin(spatial_dimension);

  Vector<Real>::const_iterator<types::RVector> normal_it =
    normal.begin(spatial_dimension);

  Vector<Real>::iterator<types::RVector>traction_end =
    tractions(el_type, ghost_type).end(spatial_dimension);

  Vector<Real>::iterator<Real>delta_max_it =
    delta_max(el_type, ghost_type).begin();

  /// compute scalars
  Real beta2 =  beta*beta;

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++delta_max_it) {

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening_it->dot(*normal_it);
    types::Vector<Real> normal_opening(spatial_dimension);
    normal_opening  = (*normal_it);
    normal_opening *= normal_opening_norm;

    types::Vector<Real> tangential_opening(spatial_dimension);
    tangential_opening  = *opening_it;
    tangential_opening -=  normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    /**
     * compute effective opening displacement
     * @f$ \delta = \sqrt{
     * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
     */
    Real delta = tangential_opening_norm;
    delta *= delta * beta2;
    delta += normal_opening_norm * normal_opening_norm;
    delta = sqrt(delta);


    /// full damage case
    if ( delta == 0) {

      /// set traction to zero
      for (UInt i = 0; i < spatial_dimension; ++i)
	(*traction_it)(i) = 0.;

    }
    /// element not fully damaged
    else {

      /**
       * Compute traction loading @f$ \mathbf{T} = 
       * e \sigma_c \frac{\delta}{\delta_c} e^{-\delta/ \delta_c}@f$
       */
      /**
       * Compute traction unloading @f$ \mathbf{T} = 
       *  \frac{t_{max}}{\delta_{max}} \delta @f$
       */
      *traction_it  = tangential_opening;
      *traction_it *= beta2;
      *traction_it += normal_opening;

      /// crack opening case
      if (delta > *delta_max_it) {
      
      Real k = exp(1)*sigma_c*exp(- delta / delta_c)/delta_c;
	*traction_it *= k;

	/// update maximum displacement
	*delta_max_it = delta;

      }

      /// unloading-reloading case
      else {

	Real k = exp(1)*sigma_c*exp(- *delta_max_it / delta_c)/delta_c;
	*traction_it *= k;
      }

    }

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::computeTangentStiffness(const ElementType & el_type,                                                          Vector<Real> & tangent_matrix,
							  const Vector<Real> & normal,
							  GhostType ghost_type) {

  
 AKANTU_DEBUG_IN();

 Vector<Real>::iterator<types::Matrix> tangent_it = tangent_matrix.begin(spatial_dimension, spatial_dimension);
 Vector<Real>::iterator<types::Matrix> tangent_end = tangent_matrix.end(spatial_dimension, spatial_dimension);
 Vector<Real>::const_iterator<types::RVector> normal_it = normal.begin(spatial_dimension);
 Vector<Real>::iterator<types::RVector> opening_it = opening(el_type, ghost_type).begin(spatial_dimension);
 Vector<Real>::iterator<types::RVector> traction_it = tractions(el_type, ghost_type).begin(spatial_dimension);
 Vector<Real>::iterator<Real>delta_max_it = delta_max(el_type, ghost_type).begin();
 Real beta2 = beta*beta;


 /**
  * compute tangent matrix  @f$ \frac{\partial \mathbf{t}} 
  * {\partial \delta} = \hat{\mathbf{t}} \otimes
  * \frac{\partial (t/\delta)}{\partial \delta}
  * \frac{\hat{\mathbf{t}}}{\delta}+ \frac{t}{\delta}  [ \beta^2 \mathbf{I} +
  * (1-\beta^2) (\mathbf{n} \otimes \mathbf{n})] @f$
  **/

 /**
  * In which @f$
  *  \frac{\partial(t/ \delta)}{\partial \delta} = 
  * \left\{\begin{array} {l l}
  *  -e  \frac{\sigma_c}{\delta_c^2  }e^{-\delta  /  \delta_c} &  \quad  if
  *  \delta \geq \delta_{max} \\
  *  0 & \quad if \delta < \delta_{max}, \delta_n > 0
  *  \end{array}\right. @f$ 
  **/

 for (; tangent_it != tangent_end; ++tangent_it, ++normal_it, ++opening_it, ++traction_it) {
   Real normal_opening_norm = opening_it->dot(*normal_it);
   types::Vector<Real> normal_opening(spatial_dimension);
   normal_opening  = (*normal_it);
   normal_opening *= normal_opening_norm;

   types::Vector<Real> tangential_opening(spatial_dimension);
   tangential_opening  = *opening_it;
   tangential_opening -=  normal_opening;

   Real tangential_opening_norm = tangential_opening.norm();

   Real delta = tangential_opening_norm;
   delta *= delta * beta2;
   delta += normal_opening_norm * normal_opening_norm;
   delta = sqrt(delta);

  
   types::RVector t_hat(tangential_opening);
   t_hat *= beta2;
   t_hat += normal_opening;
 
   types::Matrix nn(spatial_dimension, spatial_dimension);
   nn.outerProduct(*normal_it, *normal_it);

   
   types::Matrix I(spatial_dimension, spatial_dimension);
   I.eye(beta2);
   nn *= (1-beta2);
   I += nn;

   if( delta == 0.0){

     *tangent_it += I;
     *tangent_it *= exp(1)* sigma_c/delta_c;
   
   }else{


     Real traction_norm = traction_it->norm();

     I  *= traction_norm / delta;

     types::RVector t_hat_tmp (t_hat);
     Real temp_var = 0;

     if ((delta > *delta_max_it)||(abs(delta - *delta_max_it)< 1e-12)) {

       temp_var = -exp(1- delta/delta_c) * sigma_c/(delta_c * delta_c);

     }// else if((delta > 0) && (normal_opening_norm < 0)){

     //   temp_var =  exp((delta + delta_c)/delta_c) * sigma_c/(delta_c * delta_c) * (2 + delta / delta_c);

     // }
   
     temp_var /=  delta;
     t_hat_tmp *= temp_var;

     types::Matrix t_var_t(spatial_dimension, spatial_dimension);
     t_var_t.outerProduct(t_hat, t_hat_tmp);

     *tangent_it += I;
     *tangent_it += t_var_t;}
 }

 AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveExponential<spatial_dimension>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_cohesive_linear> [" << std::endl;
  stream << space << " + sigma_c      : " << sigma_c << std::endl;
  stream << space << " + beta         : " << beta << std::endl;
  stream << space << " + G_cI         : " << G_cI << std::endl;
  stream << space << " + G_cII        : " << G_cII << std::endl;
  if(this->isInit()) {
    stream << space << " + kappa      : " << kappa << std::endl;
    stream << space << " + delta_c    : " << delta_c << std::endl;
  }
  MaterialCohesive::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialCohesiveExponential);

__END_AKANTU__
