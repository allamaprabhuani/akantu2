/**
 * @file   material_cohesive_bilinear.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Wed Feb 22 16:31:20 2012
 *
 * @brief  Bilinear cohesive constitutive law
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
#include "material_cohesive_bilinear.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveBilinear<spatial_dimension>::MaterialCohesiveBilinear(SolidMechanicsModel & model, const ID & id) :
  MaterialCohesive(model,id) {
  AKANTU_DEBUG_IN();

  this->registerParam("delta_0", delta_0, 0. , _pat_parsable, "Elastic limit displacement");
  this->registerParam("beta"   , beta   , 0. , _pat_parsable, "Beta parameter"         );
  this->registerParam("G_cI"   , G_cI   , 0. , _pat_parsable, "Mode I fracture energy" );
  this->registerParam("G_cII"  , G_cII  , 0. , _pat_parsable, "Mode II fracture energy");
  this->registerParam("kappa"  , kappa  , 0. , _pat_readable, "Kappa parameter"        );
  this->registerParam("delta_c", delta_c, 0. , _pat_readable, "Critical displacement"  );
  this->registerParam("penalty", penalty, 0. , _pat_parsable, "Penalty coefficient"    );

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveBilinear<spatial_dimension>::~MaterialCohesiveBilinear() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  kappa   = G_cII / G_cI;
  delta_c = 2 * G_cI / sigma_c;

  /**
   * Recompute sigma_c as
   * @f$ {\sigma_c}_\textup{new} =
   * \frac{{\sigma_c}_\textup{old} \delta_c} {\delta_c - \delta_0} @f$
   */

  AKANTU_DEBUG_ASSERT(std::abs((this->delta_c - delta_0) / delta_0)
		      >= std::numeric_limits<Real>::epsilon(),
		      "Check your material.dat");

  this->sigma_c *= this->delta_c / (this->delta_c - delta_0);

  updateDeltaMax(_ghost);
  updateDeltaMax(_not_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::resizeCohesiveVectors() {
  MaterialCohesive::resizeCohesiveVectors();
  updateDeltaMax(_ghost);
  updateDeltaMax(_not_ghost);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::updateDeltaMax(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = this->model->getFEM("CohesiveFEM").getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);

  for(; it != last_type; ++it) {
    Vector<Real>::iterator<Real>delta_max_it =
      this->delta_max(*it, ghost_type).begin();

    Vector<Real>::iterator<Real>delta_max_end =
      this->delta_max(*it, ghost_type).end();

    for (; delta_max_it != delta_max_end; ++delta_max_it) {
      *delta_max_it = delta_0;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::computeTraction(const Vector<Real> & normal,
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

  Vector<Real>::iterator<Real>damage_it =
    damage(el_type, ghost_type).begin();

  /// compute scalars
  Real beta2_kappa2 = beta*beta/kappa/kappa;
  Real beta2_kappa  = beta*beta/kappa;

  Real epsilon = std::numeric_limits<Real>::epsilon();

  Real * memory_space = new Real[2*spatial_dimension];
  types::Vector<Real> normal_opening(memory_space, spatial_dimension);
  types::Vector<Real> tangential_opening(memory_space + spatial_dimension,
					 spatial_dimension);

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++delta_max_it, ++damage_it) {

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
    if (delta >= delta_c || delta <= epsilon) {
      /// set traction to zero
      (*traction_it).clear();

      if (normal_opening_norm < 0) {
      	normal_opening *= penalty;
      	*traction_it += normal_opening;
      }
      else {
	*damage_it = delta >= delta_c;
	*delta_max_it = *damage_it * delta_c;
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
      *damage_it = *delta_max_it / delta_c;

      Real k = sigma_c / *delta_max_it * (1. - *damage_it);
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

INSTANSIATE_MATERIAL(MaterialCohesiveBilinear);


__END_AKANTU__
