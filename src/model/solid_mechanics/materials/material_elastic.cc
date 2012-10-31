/**
 * @file   material_elastic.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the elastic material
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
#include "material_elastic.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialElastic<spatial_dimension>::MaterialElastic(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("E"           ,E           , 0.   , ParamAccessType(_pat_parsable | _pat_modifiable), "Young's modulus"        );
  this->registerParam("nu"          ,nu          , 0.5  , ParamAccessType(_pat_parsable | _pat_modifiable), "Poisson's ratio"        );
  this->registerParam("Plane_Stress",plane_stress, false, ParamAccessType(_pat_parsable | _pat_readable)  , "Is plane stress"        );
  this->registerParam("lambda"      ,lambda             , _pat_readable, "First Lamé coefficient" );
  this->registerParam("mu"          ,mu                 , _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa"        ,kpa                , _pat_readable, "Bulk coefficient"       );

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();
  if (spatial_dimension == 1) nu = 0.;
  updateInternalParameters();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::updateInternalParameters() {
  lambda   = nu * E / ((1 + nu) * (1 - 2*nu));
  mu       = E / (2 * (1 + nu));

  if(plane_stress) lambda = nu * E / ((1 + nu)*(1 - nu));

  kpa      = lambda + 2./3. * mu;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  computeStressOnQuad(grad_u, sigma);
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
							      Vector<Real> & tangent_matrix,
							      __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  computeTangentModuliOnQuad(tangent);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElastic<spatial_dimension>::getPushWaveSpeed() const {
  return sqrt((lambda + 2*mu)/rho);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElastic<spatial_dimension>::getShearWaveSpeed() const {
  return sqrt(mu/rho);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialElastic);

__END_AKANTU__
