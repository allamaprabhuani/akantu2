/**
 * @file   material_cohesive_linear.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Feb 20 12:14:16 2012
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading
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
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialCohesiveLinear::MaterialCohesiveLinear(Model & model, const ID & id) :
  MaterialCohesive(model,id),
  delta_max("delta max",id) {
  AKANTU_DEBUG_IN();

  beta      = 0;
  G_cI      = 0;
  G_cII     = 0;

  initInternalVector(delta_max, 1, _ek_cohesive);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialCohesiveLinear::~MaterialCohesiveLinear() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesiveLinear::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  kappa   = G_cII / G_cI;
  delta_c = 2 * G_cI / sigma_c;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesiveLinear::resizeCohesiveVectors() {
  MaterialCohesive::resizeCohesiveVectors();
  resizeInternalVector(delta_max, _ek_cohesive);
}

/* -------------------------------------------------------------------------- */
bool MaterialCohesiveLinear::setParam(const std::string & key, 
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
void MaterialCohesiveLinear::computeTraction(const Vector<Real> & normal,
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
  Real beta2_kappa2 = beta*beta/kappa/kappa;
  Real beta2_kappa  = beta*beta/kappa;


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
    delta *= delta * beta2_kappa2;
    delta += normal_opening_norm * normal_opening_norm;
    delta = sqrt(delta);


    /// full damage case
    if (delta >= delta_c || delta == 0) {

      /// set traction to zero
      (*traction_it).clear();

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
      *traction_it += normal_opening;

      /// update maximum displacement
      *delta_max_it = std::max(*delta_max_it, delta);

      Real k = sigma_c / *delta_max_it * (1. - *delta_max_it / delta_c);
      *traction_it *= k;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesiveLinear::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_cohesive_linear> [" << std::endl;
  stream << space << " + sigma_c      : " << sigma_c << std::endl;
  stream << space << " + beta         : " << beta << std::endl;
  stream << space << " + G_cI         : " << G_cI << std::endl;
  stream << space << " + G_cII        : " << G_cII << std::endl;
  if(is_init) {
    stream << space << " + kappa      : " << kappa << std::endl;
    stream << space << " + delta_c    : " << delta_c << std::endl;
  }
  MaterialCohesive::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */


__END_AKANTU__
