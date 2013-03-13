/**
 * @file   material_cohesive_linear_exponential_extrinsic.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Thu May 24 10:46:59 2012
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
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
#include "material_cohesive_linear_exponential_extrinsic.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveLinearExponentialExtrinsic<spatial_dimension>::MaterialCohesiveLinearExponentialExtrinsic(SolidMechanicsModel & model, const ID & id) :
  MaterialCohesive(model,id),
  sigma_c_eff("sigma_c_eff",id),
  sigma_actual("sigma actual",id) {
  AKANTU_DEBUG_IN();

  sigma_max = 0;
  delta_0   = 0;
  beta      = 0;
  G_cI      = 0;
  G_cII     = 0;

  initInternalArray(sigma_c_eff, 1, _ek_cohesive);
  initInternalArray(sigma_actual, 1, _ek_cohesive);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveLinearExponentialExtrinsic<spatial_dimension>::~MaterialCohesiveLinearExponentialExtrinsic() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearExponentialExtrinsic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  kappa   = G_cII / G_cI;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearExponentialExtrinsic<spatial_dimension>::resizeCohesiveArrays() {
  MaterialCohesive::resizeCohesiveArrays();

  resizeInternalArray(sigma_c_eff, _ek_cohesive);
  resizeInternalArray(sigma_actual, _ek_cohesive);

  FEM & fem_cohesive = model->getFEM("CohesiveFEM");
  const Mesh & mesh = fem_cohesive.getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);
  for(; it != last_type; ++it) {

    const Array<UInt> & elem_filter = element_filter(*it, _not_ghost);
    UInt nb_element = elem_filter.getSize();

    if (nb_element == 0) continue;

    UInt nb_quadrature_points = fem_cohesive.getNbQuadraturePoints(*it, _not_ghost);
    UInt nb_element_old = nb_element - sigma_insertion.getSize() / nb_quadrature_points;

    Array<Real> & sigma_c_eff_vec = sigma_c_eff(*it, _not_ghost);
    Array<Real> & sigma_actual_vec = sigma_actual(*it, _not_ghost);

    for (UInt el = nb_element_old; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quadrature_points; ++q) {
	Real new_sigma = sigma_insertion((el - nb_element_old)*nb_quadrature_points + q);
	sigma_c_eff_vec(el * nb_quadrature_points + q) = new_sigma;
	sigma_actual_vec(el * nb_quadrature_points + q) = new_sigma;
      }
    }
  }

}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialCohesiveLinearExponentialExtrinsic<spatial_dimension>::computeEffectiveNorm(const Matrix<Real> & stress, const Vector<Real> & normal, const Vector<Real> & tangent) {
  AKANTU_DEBUG_IN();

  Real normal_contrib, tangent_contrib;
  Vector<Real> normal_stress(spatial_dimension);
  Vector<Real> tangential_stress(spatial_dimension);

  normal_stress.mul<false>(stress, normal);
  tangential_stress.mul<false>(stress, tangent);

  normal_contrib = normal_stress.dot(normal);
  tangent_contrib = tangential_stress.dot(tangent);

  if (normal_contrib < 0) normal_contrib = 0;

  AKANTU_DEBUG_OUT();

  Real epsilon = std::numeric_limits<Real>::epsilon();

  if (std::abs(beta) < epsilon)
    return normal_contrib;
  else
    return std::sqrt(normal_contrib*normal_contrib
		     + tangent_contrib*tangent_contrib/beta/beta);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearExponentialExtrinsic<spatial_dimension>::computeTraction(const Array<Real> & normal,
								 ElementType el_type,
								 GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  Array<Real>::iterator< Vector<Real> > traction_it =
    tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> > opening_it =
    opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::const_iterator< Vector<Real> > normal_it =
    normal.begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> >traction_end =
    tractions(el_type, ghost_type).end(spatial_dimension);

  Array<Real>::iterator<Real>sigma_c_it =
    sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>delta_max_it =
    delta_max(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>sigma_actual_it =
    sigma_actual(el_type, ghost_type).begin();

  /// compute scalars
  Real beta2_kappa2 = beta*beta/kappa/kappa;
  Real beta2_kappa  = beta*beta/kappa;

  Real epsilon = std::numeric_limits<Real>::epsilon();

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++sigma_c_it,
	 ++delta_max_it, ++sigma_actual_it) {

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening_it->dot(*normal_it);
    Vector<Real> normal_opening(spatial_dimension);
    normal_opening  = (*normal_it);
    normal_opening *= normal_opening_norm;

    Vector<Real> tangential_opening(spatial_dimension);
    tangential_opening  = *opening_it;
    tangential_opening -=  normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    /**
     * compute effective opening displacement
     * @f$ \delta = \sqrt{
     * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
     */
    Real delta = tangential_opening_norm;
    delta *= delta * beta2_kappa2;
    delta += normal_opening_norm * normal_opening_norm;
    delta = sqrt(delta);

    *traction_it  = tangential_opening;
    *traction_it *= beta2_kappa;
    *traction_it += normal_opening;
    *traction_it /= delta;

    /// crack opening case
    if (std::abs(delta) <= std::abs(delta) * epsilon || normal_opening_norm < 0) {
      (*traction_it).clear();
    }
    else if (delta > *delta_max_it) {

      Real sigma;

      if (delta < delta_0) {
	sigma = *sigma_c_it + (sigma_max - *sigma_c_it) / delta_0 * delta;
      }
      else {
	Real z = -21.61440010 + 0.04166666667 * sqrt(2.69097e5 + 4.8e5 * delta);
	z = std::min(z, z_max);
	sigma = exp(-1*gamma*z) * (sigma_max - sigma_max / z_max * z);
      }

      *traction_it *= sigma;

      /// update maximum displacement and stress
      *delta_max_it = delta;
      *sigma_actual_it = sigma;
    }
    /// unloading-reloading case
    else {
      *traction_it *= *sigma_actual_it / *delta_max_it * delta;
    }

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialCohesiveLinearExponentialExtrinsic);

__END_AKANTU__
