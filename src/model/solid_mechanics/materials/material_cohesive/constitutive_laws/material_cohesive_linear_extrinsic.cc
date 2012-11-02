/**
 * @file   material_cohesive_linear.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Feb 20 12:14:16 2012
 *
 * @brief Linear irreversible cohesive law of mixed mode loading with
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
#include "material_cohesive_linear_extrinsic.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveLinearExtrinsic<spatial_dimension>::MaterialCohesiveLinearExtrinsic(SolidMechanicsModel & model,
										    const ID & id) :
  MaterialCohesive(model,id),
  sigma_c_eff("sigma_c_eff",id),
  delta_c("delta_c",id),
  normal_stress(spatial_dimension),
  tangential_stress(spatial_dimension) {
  AKANTU_DEBUG_IN();

  this->registerParam("beta"   , beta   , 0. , _pat_parsable, "Beta parameter"         );
  this->registerParam("G_cI"   , G_cI   , 0. , _pat_parsable, "Mode I fracture energy" );
  this->registerParam("G_cII"  , G_cII  , 0. , _pat_parsable, "Mode II fracture energy");
  this->registerParam("penalty", penalty, 0. , _pat_parsable, "Penalty coefficient"    );
  this->registerParam("kappa"  , kappa  , 1. , _pat_readable, "Kappa parameter"        );

  initInternalVector(sigma_c_eff, 1, false, _ek_cohesive);
  initInternalVector(    delta_c, 1, false, _ek_cohesive);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveLinearExtrinsic<spatial_dimension>::~MaterialCohesiveLinearExtrinsic() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearExtrinsic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  if (G_cII != 0)
    kappa = G_cII / G_cI;

  /// compute scalars
  beta2_kappa2 = beta * beta/kappa/kappa;
  beta2_kappa  = beta * beta/kappa;

  if (beta == 0)
    beta2_inv = 0;
  else
    beta2_inv = 1./beta/beta;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearExtrinsic<spatial_dimension>::resizeCohesiveVectors() {
  MaterialCohesive::resizeCohesiveVectors();

   resizeInternalVector(sigma_c_eff, _ek_cohesive);
   resizeInternalVector(delta_c, _ek_cohesive);

  FEM & fem_cohesive = model->getFEM("CohesiveFEM");
  const Mesh & mesh = fem_cohesive.getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);
  for(; it != last_type; ++it) {

    const Vector<UInt> & elem_filter = element_filter(*it, _not_ghost);
    UInt nb_element = elem_filter.getSize();

    if (nb_element == 0) continue;

    UInt nb_quadrature_points = fem_cohesive.getNbQuadraturePoints(*it, _not_ghost);
    UInt nb_element_old = nb_element - sigma_insertion.getSize() / nb_quadrature_points;

    Vector<Real> & sigma_c_eff_vec = sigma_c_eff(*it, _not_ghost);
    Vector<Real> & delta_c_vec = delta_c(*it, _not_ghost);

    for (UInt el = nb_element_old; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quadrature_points; ++q) {
	Real new_sigma = sigma_insertion((el - nb_element_old) * nb_quadrature_points+q);
	Real new_delta = 2 * G_cI / new_sigma;
	sigma_c_eff_vec(el * nb_quadrature_points + q) = new_sigma;
	delta_c_vec(el * nb_quadrature_points + q) = new_delta;
      }
    }
  }

}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialCohesiveLinearExtrinsic<spatial_dimension>::computeEffectiveNorm(const types::Matrix & stress,
										     const types::RVector & normal,
										     const types::RVector & tangent) {
  AKANTU_DEBUG_IN();

  normal_stress.mul<false>(stress, normal);
  tangential_stress.mul<false>(stress, tangent);

  Real normal_contrib = normal_stress.dot(normal);
  Real tangent_contrib = tangential_stress.dot(tangent);

  normal_contrib = std::max(0., normal_contrib);

  AKANTU_DEBUG_OUT();

  return std::sqrt(normal_contrib*normal_contrib
		   + tangent_contrib*tangent_contrib*beta2_inv);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearExtrinsic<spatial_dimension>::computeStressNorms(const Vector<Real> & facet_stress,
									    Vector<Real> & stress_check) {
  AKANTU_DEBUG_IN();

  sigma_insertion.resize(0);

  Vector<bool> & facets_check = model->getFacetsCheck();
  ElementType type_facet = model->getFacetType();

  UInt nb_quad_facet = model->getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
  UInt nb_facet = facets_check.getSize();

  const Vector<Real> & tangents = model->getTangents();
  const Vector<Real> & normals
    = model->getFEM("FacetsFEM").getNormalsOnQuadPoints(type_facet);

  Vector<Real>::iterator<types::RVector> stress_check_it =
    stress_check.begin(nb_quad_facet);

  Vector<Real>::const_iterator<types::RVector> normal_it =
    normals.begin(spatial_dimension);

  Vector<Real>::const_iterator<types::RVector> tangent_it =
    tangents.begin(spatial_dimension);

  Vector<Real>::const_iterator<types::Matrix> facet_stress_it =
    facet_stress.begin(spatial_dimension, spatial_dimension);

  for (UInt f = 0; f < nb_facet; ++f, ++stress_check_it) {
    for (UInt q = 0; q < nb_quad_facet; ++q, ++normal_it, ++tangent_it) {
      for (UInt e = 0; e < 2; ++e, ++facet_stress_it) {

	if (facets_check(f) == true) {
	  Real effective_norm =
	    computeEffectiveNorm(*facet_stress_it, *normal_it, *tangent_it);

	  (*stress_check_it)(q) =
	    std::max((*stress_check_it)(q), effective_norm);
	}

      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinearExtrinsic<spatial_dimension>::computeTraction(const Vector<Real> & normal,
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

  Vector<Real>::iterator<Real>sigma_c_it =
    sigma_c_eff(el_type, ghost_type).begin();

  Vector<Real>::iterator<Real>delta_max_it =
    delta_max(el_type, ghost_type).begin();

  Vector<Real>::iterator<Real>delta_c_it =
    delta_c(el_type, ghost_type).begin();

  Vector<Real>::iterator<Real>damage_it =
    damage(el_type, ghost_type).begin();

  Real epsilon = std::numeric_limits<Real>::epsilon();

  Real * memory_space = new Real[2*spatial_dimension];
  types::Vector<Real> normal_opening(memory_space, spatial_dimension);
  types::Vector<Real> tangential_opening(memory_space + spatial_dimension,
					 spatial_dimension);

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++sigma_c_it,
	 ++delta_max_it, ++delta_c_it, ++damage_it) {

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening_it->dot(*normal_it);
    normal_opening  = (*normal_it);
    normal_opening *= normal_opening_norm;

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

    /// don't consider penetration contribution
    if (normal_opening_norm > 0)
      delta += normal_opening_norm * normal_opening_norm;

    delta = sqrt(delta);

    /// full damage case or zero displacement case
    if (delta >= *delta_c_it || delta <= epsilon) {

      /// set traction to zero
      (*traction_it).clear();

      if (normal_opening_norm < 0) {
      	normal_opening *= penalty;
      	*traction_it += normal_opening;
      }
      else {
	*damage_it = delta >= *delta_c_it;
	*delta_max_it = *damage_it * (*delta_c_it);
      }
    }
    /// element not fully damaged
    else {

      /**
       * Compute traction @f$ \mathbf{T} = \left(
       * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
       * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
       * \frac{\delta}{\delta_c} \right)@f$
       */

      *traction_it  = tangential_opening;
      *traction_it *= beta2_kappa;

      /// update maximum displacement
      *delta_max_it = std::max(*delta_max_it, delta);
      *damage_it = *delta_max_it / *delta_c_it;

      Real k = *sigma_c_it / *delta_max_it * (1. - *damage_it);
      *traction_it *= k;

      /// use penalty coefficient in case of penetration
      if (normal_opening_norm < 0) k = penalty;

      normal_opening *= k;
      *traction_it += normal_opening;
    }
  }

  delete [] memory_space;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialCohesiveLinearExtrinsic);

__END_AKANTU__
