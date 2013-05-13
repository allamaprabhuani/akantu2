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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveBilinear<spatial_dimension>::MaterialCohesiveBilinear(SolidMechanicsModel & model, const ID & id) :
  MaterialCohesiveLinear<spatial_dimension>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("delta_0", delta_0, 0. , _pat_parsable, "Elastic limit displacement");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  /**
   * Recompute sigma_c as
   * @f$ {\sigma_c}_\textup{new} =
   * \frac{{\sigma_c}_\textup{old} \delta_c} {\delta_c - \delta_0} @f$
   */

  current_delta_c = 2 * this->G_cI / this->sigma_c;

  AKANTU_DEBUG_ASSERT(current_delta_c > delta_0, "Check your material file");

  this->sigma_c *= current_delta_c / (current_delta_c - delta_0);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::resizeCohesiveArrays() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::resizeCohesiveArrays();
  this->resizeInternalArray(this->sigma_c_eff, _ek_cohesive);
  this->resizeInternalArray(this->delta_c, _ek_cohesive);

  const Mesh & mesh = this->model->getFEM("CohesiveFEM").getMesh();

  GhostType gt = _not_ghost;
  ElementKind element_kind = _ek_cohesive;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, element_kind);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, element_kind);
  for(; it != end; ++it) {
    this->sigma_c_eff(*it, gt).set(this->sigma_c);
    this->delta_c(*it, gt).set(current_delta_c);
    this->delta_max(*it, gt).set(delta_0);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialCohesiveBilinear);


__END_AKANTU__
