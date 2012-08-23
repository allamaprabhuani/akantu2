/**
 * @file   material_marigo_non_local_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue May  1 14:12:17 2012
 *
 * @brief  
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

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
MaterialMarigoNonLocal<spatial_dimension, WeigthFunction>::MaterialMarigoNonLocal(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id),
  MaterialMarigoNonLocalParent(model, id),
  Y("Y", id), Ynl("Y non local", id) {
  AKANTU_DEBUG_IN();
  this->is_non_local = true;
  this->initInternalVector(this->Y, 1);
  this->initInternalVector(this->Ynl, 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialMarigoNonLocal<spatial_dimension, WeigthFunction>::initMaterial() {
  AKANTU_DEBUG_IN();

  this->resizeInternalVector(this->Y);
  this->registerNonLocalVariable(Y, Ynl, 1);

  MaterialMarigoNonLocalParent::initMaterial();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialMarigoNonLocal<spatial_dimension, WeigthFunction>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Yt  = this->Y(el_type, ghost_type).storage();
  Real * Ydq = this->Yd_rand(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  MaterialMarigo<spatial_dimension>::computeStressOnQuad(grad_u, sigma,
							 *dam, *Yt, *Ydq);
  ++dam;
  ++Yt;
  ++Ydq;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialMarigoNonLocal<spatial_dimension, WeigthFunction>::computeNonLocalStress(ElementType type,
										      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam  = this->damage(type, ghost_type).storage();
  Real * Ydq  = this->Yd_rand(type, ghost_type).storage();
  Real * Ynlt = this->Ynl(type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(type, ghost_type);
  this->computeDamageAndStressOnQuad(sigma, *dam, *Ynlt, *Ydq);

  ++dam;
  ++Ynlt;
  ++Ydq;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}
