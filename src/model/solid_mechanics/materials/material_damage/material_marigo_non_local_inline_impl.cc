/**
 * @file   material_marigo_non_local_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  Marigo non-local inline function implementation
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
__END_AKANTU__

#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#include <string>
#endif

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialMarigoNonLocal<spatial_dimension>::MaterialMarigoNonLocal(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id),
  MaterialMarigoNonLocalParent(model, id),
  Y("Y", *this), Ynl("Y non local", *this) {
  AKANTU_DEBUG_IN();
  this->is_non_local = true;
  this->Y.initialize(1);
  this->Ynl.initialize(1);
  this->model->getNonLocalManager().registerNonLocalVariable(this->Y.getName(), Ynl.getName(), 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMarigoNonLocal<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialMarigoNonLocalParent::initMaterial();
  this->model->getNonLocalManager().nonLocalVariableToNeighborhood(Ynl.getName(), this->name);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMarigoNonLocal<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Yt  = this->Y(el_type, ghost_type).storage();
  Real * Ydq = this->Yd(el_type, ghost_type).storage();

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
template<UInt spatial_dimension>
void MaterialMarigoNonLocal<spatial_dimension>::computeNonLocalStress(ElementType type,
										      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam  = this->damage(type, ghost_type).storage();
  Real * Ydq  = this->Yd(type, ghost_type).storage();
  Real * Ynlt = this->Ynl(type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(type, ghost_type);
  this->computeDamageAndStressOnQuad(sigma, *dam, *Ynlt, *Ydq);

  ++dam;
  ++Ynlt;
  ++Ydq;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

