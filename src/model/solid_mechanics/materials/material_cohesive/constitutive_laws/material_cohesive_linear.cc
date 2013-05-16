/**
 * @file   material_cohesive_linear_extrinsic.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Tue May 08 13:01:18 2012
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
#include "material_cohesive_linear.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveLinear<spatial_dimension>::MaterialCohesiveLinear(SolidMechanicsModel & model,
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

  AKANTU_DEBUG_OUT();
  }

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::initMaterial() {
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

  initInternalArray(sigma_c_eff, 1, false, _ek_cohesive);
  initInternalArray(    delta_c, 1, false, _ek_cohesive);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::checkInsertion(const ByElementTypeReal & facet_stress,
							       const Mesh & mesh_facets,
							       ByElementTypeArray<bool> & facet_insertion) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1);

  for (; it != last; ++it) {

    ElementType type_facet = *it;
    ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);
    Array<bool> & facets_check = model->getFacetsCheck();
    Array<bool> & f_insertion = facet_insertion(type_facet);
    Array<UInt> & f_filter = facet_filter(type_facet);
    Array<Real> & sig_c_eff = sigma_c_eff(type_cohesive);
    Array<Real> & del_c = delta_c(type_cohesive);
    const Array<Real> & f_stress = facet_stress(type_facet);
    const Array<Real> & sigma_lim = sigma_limit(type_facet);

    UInt nb_quad_facet = model->getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
    UInt nb_facet = f_filter.getSize();
    if (nb_facet == 0) continue;

    Array<Real> stress_check(nb_facet, nb_quad_facet);
    stress_check.clear();

    computeStressNorms(f_stress, stress_check, type_facet);

    Array<Real>::iterator< Vector<Real> > stress_check_it =
      stress_check.begin(nb_quad_facet);

    Array<Real>::const_iterator<Real> sigma_lim_it = sigma_lim.begin();

    for (UInt f = 0; f < nb_facet;
	 ++f, ++stress_check_it, ++sigma_lim_it) {
      UInt facet = f_filter(f);
      if (facets_check(facet) == true) {

	for (UInt q = 0; q < nb_quad_facet; ++q) {
	  if ((*stress_check_it)(q) > (*sigma_lim_it)) {
	    f_insertion(facet) = true;
	    for (UInt qs = 0; qs < nb_quad_facet; ++qs) {
	      Real new_sigma = (*stress_check_it)(qs);
	      Real new_delta = 2 * G_cI / new_sigma;
	      sig_c_eff.push_back(new_sigma);
	      del_c.push_back(new_delta);
	    }
	    facets_check(facet) = false;
	    break;
	  }
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialCohesiveLinear<spatial_dimension>::computeEffectiveNorm(const Matrix<Real> & stress,
									    const Vector<Real> & normal,
									    const Vector<Real> & tangent) {
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
void MaterialCohesiveLinear<spatial_dimension>::computeStressNorms(const Array<Real> & facet_stress,
								   Array<Real> & stress_check,
								   ElementType type_facet) {
  AKANTU_DEBUG_IN();

  Array<bool> & facets_check = model->getFacetsCheck();
  Array<UInt> & f_filter = facet_filter(type_facet);

  UInt nb_quad_facet = model->getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);

  const Array<Real> & tangents = model->getTangents(type_facet);
  const Array<Real> & normals
    = model->getFEM("FacetsFEM").getNormalsOnQuadPoints(type_facet);

  UInt nb_facet = normals.getSize();

  Array<Real>::iterator< Vector<Real> > stress_check_it =
    stress_check.begin(nb_quad_facet);

  Array<Real>::const_iterator< Vector<Real> > normal_it =
    normals.begin(spatial_dimension);

  Array<Real>::const_iterator< Vector<Real> > tangent_it =
    tangents.begin(spatial_dimension);

  Array<Real>::const_iterator< Matrix<Real> > facet_stress_it =
    facet_stress.begin(spatial_dimension, spatial_dimension);

  UInt facet_index = 0;
  UInt facet = f_filter(facet_index);
  UInt nq2 = nb_quad_facet * 2;

  for (UInt f = 0; f < nb_facet; ++f) {

    if (f == facet) {
      for (UInt q = 0; q < nb_quad_facet; ++q, ++normal_it, ++tangent_it) {
	for (UInt e = 0; e < 2; ++e, ++facet_stress_it) {

	  if (facets_check(facet) == true) {
	    Real effective_norm =
	      computeEffectiveNorm(*facet_stress_it, *normal_it, *tangent_it);

	    (*stress_check_it)(q) =
	      std::max((*stress_check_it)(q), effective_norm);
	  }

	}
      }
      ++facet_index;
      if (facet_index == f_filter.getSize()) break;
      ++stress_check_it;
      facet = f_filter(facet_index);
    }
    else {
      normal_it += nb_quad_facet;
      tangent_it += nb_quad_facet;
      facet_stress_it += nq2;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::computeTraction(const Array<Real> & normal,
								ElementType el_type,
								GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  Array<Real>::iterator< Vector<Real> > traction_it =
    tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> > opening_it =
    opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> > contact_traction_it =
    contact_tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::const_iterator< Vector<Real> > normal_it =
    normal.begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> >traction_end =
    tractions(el_type, ghost_type).end(spatial_dimension);

  Array<Real>::iterator<Real>sigma_c_it =
    sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>delta_max_it =
    delta_max(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>delta_c_it =
    delta_c(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>damage_it =
    damage(el_type, ghost_type).begin();

  Real * memory_space = new Real[2*spatial_dimension];
  Vector<Real> normal_opening(memory_space, spatial_dimension);
  Vector<Real> tangential_opening(memory_space + spatial_dimension,
				  spatial_dimension);

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++sigma_c_it,
	 ++delta_max_it, ++delta_c_it, ++damage_it, ++contact_traction_it) {

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
    Real delta = tangential_opening_norm * tangential_opening_norm;
    delta *= beta2_kappa2;

    /// don't consider penetration contribution
    if (normal_opening_norm > 0)
      delta += normal_opening_norm * normal_opening_norm;

    delta = std::sqrt(delta);

    /// use penalty coefficient in case of penetration
    contact_traction_it->clear();
    if (normal_opening_norm < 0) {
      *contact_traction_it = normal_opening;
      *contact_traction_it *= penalty;
      normal_opening.set(0.);
    }

    /// update maximum displacement and damage
    *delta_max_it = std::max(*delta_max_it, delta);
    *damage_it = std::min(*delta_max_it / *delta_c_it, 1.);

    /**
     * Compute traction @f$ \mathbf{T} = \left(
     * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
     * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
     * \frac{\delta}{\delta_c} \right)@f$
     */

    if (Math::are_float_equal(*delta_max_it, 0.) ||
	Math::are_float_equal(*damage_it, 1.))
      traction_it->clear();
    else {
      *traction_it  = tangential_opening;
      *traction_it *= beta2_kappa;
      *traction_it += normal_opening;

      *traction_it *= *sigma_c_it / *delta_max_it * (1. - *damage_it);
    }
  }

  delete [] memory_space;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialCohesiveLinear);

__END_AKANTU__
