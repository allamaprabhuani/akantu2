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
  insertion_stress("insertion_stress", id) {
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

  initInternalArray(      sigma_c_eff,                 1, false, _ek_cohesive);
  initInternalArray(          delta_c,                 1, false, _ek_cohesive);
  initInternalArray( insertion_stress, spatial_dimension, false, _ek_cohesive);

  /// keep tolerance small enough for the constitutive law
  if (sigma_c != 0) {
    Real delta_limit = 2 * G_cI / sigma_c / 20.;
    if (Math::getTolerance() > delta_limit)
      Math::setTolerance(delta_limit);
  }

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
    Array<Real> & ins_stress = insertion_stress(type_cohesive);
    const Array<Real> & f_stress = facet_stress(type_facet);
    const Array<Real> & sigma_lim = sigma_limit(type_facet);

    UInt nb_quad_facet = model->getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
    UInt nb_facet = f_filter.getSize();
    if (nb_facet == 0) continue;

    UInt tot_nb_quad = nb_facet * nb_quad_facet;

    Array<Real> stress_check(tot_nb_quad);
    Array<Real> normal_stress(tot_nb_quad, spatial_dimension);

    computeStressNorms(f_stress, stress_check, normal_stress, type_facet);

    Array<Real>::iterator<Vector<Real> > stress_check_it =
      stress_check.begin_reinterpret(nb_quad_facet, nb_facet);
    Array<Real>::iterator<Matrix<Real> > normal_stress_it =
      normal_stress.begin_reinterpret(nb_quad_facet, spatial_dimension, nb_facet);

    Array<Real>::const_iterator<Real> sigma_lim_it = sigma_lim.begin();

    for (UInt f = 0; f < nb_facet; ++f, ++sigma_lim_it, ++stress_check_it,
	   ++normal_stress_it) {
      UInt facet = f_filter(f);
      if (facets_check(facet) == true) {

	for (UInt q = 0; q < nb_quad_facet; ++q) {
	  if ((*stress_check_it)(q) > (*sigma_lim_it)) {
	    f_insertion(facet) = true;
	    for (UInt qs = 0; qs < nb_quad_facet; ++qs) {
	      Real new_sigma = (*stress_check_it)(qs);
	      Real new_delta = 2 * G_cI / new_sigma;

	      Vector<Real> ins_s(normal_stress_it->storage() + qs * spatial_dimension,
				 spatial_dimension);
	      ins_s *= -1.;

	      sig_c_eff.push_back(new_sigma);
	      del_c.push_back(new_delta);
	      ins_stress.push_back(ins_s);
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
inline void MaterialCohesiveLinear<spatial_dimension>::computeEffectiveNorm(const Matrix<Real> & stress,
									    const Vector<Real> & normal,
									    const Vector<Real> & tangent,
									    Vector<Real> & normal_stress,
									    Real & effective_norm) {
  AKANTU_DEBUG_IN();

  normal_stress.mul<false>(stress, normal);

  Real normal_contrib = normal_stress.dot(normal);
  Real tangent_contrib = normal_stress.dot(tangent);

  normal_contrib = std::max(0., normal_contrib);

  effective_norm = std::sqrt(normal_contrib * normal_contrib
			     + tangent_contrib * tangent_contrib * beta2_inv);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::computeStressNorms(const Array<Real> & facet_stress,
								   Array<Real> & stress_check,
								   Array<Real> & normal_stress,
								   ElementType type_facet) {
  AKANTU_DEBUG_IN();

  Array<bool> & facets_check = model->getFacetsCheck();
  Array<UInt> & f_filter = facet_filter(type_facet);

  UInt nb_quad_facet = model->getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);

  const Array<Real> & tangents = model->getTangents(type_facet);
  const Array<Real> & normals
    = model->getFEM("FacetsFEM").getNormalsOnQuadPoints(type_facet);

  UInt nb_facet = normals.getSize();

  Array<Real>::iterator<Real > stress_check_it = stress_check.begin();

  Array<Real>::const_iterator< Vector<Real> > normal_it =
    normals.begin(spatial_dimension);

  Array<Real>::const_iterator< Vector<Real> > tangent_it =
    tangents.begin(spatial_dimension);

  Array<Real>::const_iterator< Matrix<Real> > facet_stress_it =
    facet_stress.begin(spatial_dimension, spatial_dimension);

  Array<Real>::iterator<Vector<Real> > normal_stress_it =
    normal_stress.begin(spatial_dimension);

  UInt facet_index = 0;
  UInt facet = f_filter(facet_index);
  UInt nq2 = nb_quad_facet * 2;
  Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);

  for (UInt f = 0; f < nb_facet; ++f) {

    if (f == facet) {
      for (UInt q = 0; q < nb_quad_facet; ++q, ++normal_it, ++tangent_it,
	     ++stress_check_it, ++normal_stress_it) {

	if (facets_check(facet) == true) {

	  /// compute average stress
	  stress_tmp.clear();

	  for (UInt p = 0; p < 2; ++p, ++facet_stress_it)
	    stress_tmp += *facet_stress_it;

	  stress_tmp /= 2.;

	  /// compute normal, tangential and effective stress
	  computeEffectiveNorm(stress_tmp, *normal_it, *tangent_it,
			       *normal_stress_it, *stress_check_it);
	}
	else
	  facet_stress_it += 2;

      }
      ++facet_index;
      if (facet_index == f_filter.getSize()) break;
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

  Array<Real>::iterator< Vector<Real> > contact_opening_it =
    contact_opening(el_type, ghost_type).begin(spatial_dimension);

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

  Array<Real>::iterator<Vector<Real> > insertion_stress_it =
    insertion_stress(el_type, ghost_type).begin(spatial_dimension);


  Real * memory_space = new Real[2*spatial_dimension];
  Vector<Real> normal_opening(memory_space, spatial_dimension);
  Vector<Real> tangential_opening(memory_space + spatial_dimension,
				  spatial_dimension);

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++sigma_c_it,
	 ++delta_max_it, ++delta_c_it, ++damage_it, ++contact_traction_it,
	 ++insertion_stress_it, ++contact_opening_it) {

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

    /// don't consider penetration contribution for delta
    if (normal_opening_norm > 0) {
      delta += normal_opening_norm * normal_opening_norm;
      contact_traction_it->clear();
      contact_opening_it->clear();
    }
    else {
      /// use penalty coefficient in case of penetration
      *contact_traction_it = normal_opening;
      *contact_traction_it *= penalty;
      *contact_opening_it = normal_opening;
      *opening_it = tangential_opening;
      normal_opening.clear();
    }

    delta = std::sqrt(delta);

    /// update maximum displacement and damage
    *delta_max_it = std::max(*delta_max_it, delta);
    *damage_it = std::min(*delta_max_it / *delta_c_it, 1.);

    /**
     * Compute traction @f$ \mathbf{T} = \left(
     * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
     * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
     * \frac{\delta}{\delta_c} \right)@f$
     */

    if (Math::are_float_equal(*damage_it, 1.))
      traction_it->clear();
    else if (Math::are_float_equal(*damage_it, 0.)) {
      if (Math::are_float_equal(normal_opening_norm, 0.))
	*traction_it = *insertion_stress_it;
      else
	traction_it->clear();
    }
    else {
      *traction_it  = tangential_opening;
      *traction_it *= beta2_kappa;
      *traction_it += normal_opening;

      AKANTU_DEBUG_ASSERT(*delta_max_it != 0.,
			  "Division by zero, tolerance might be too low");

      *traction_it *= *sigma_c_it / *delta_max_it * (1. - *damage_it);
    }
  }

  delete [] memory_space;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialCohesiveLinear);

__END_AKANTU__
