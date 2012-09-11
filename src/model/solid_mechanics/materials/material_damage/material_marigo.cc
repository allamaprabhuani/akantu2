/**
 * @file   material_marigo.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the marigo material
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
#include "material_marigo.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialMarigo<spatial_dimension>::MaterialMarigo(SolidMechanicsModel & model,
						  const ID & id)  :
  Material(model, id),
  MaterialDamage<spatial_dimension>(model, id),
  Yd_rand("Yd_rand",id), damage_in_y(false), yc_limit(false) {
  AKANTU_DEBUG_IN();

  this->registerParam("Yd",                       Yd,   50., ParamAccessType(_pat_parsable | _pat_modifiable));
  this->registerParam("Sd",                       Sd, 5000., ParamAccessType(_pat_parsable | _pat_modifiable));
  this->registerParam("Yd_randomness", Yd_randomness,    0., _pat_parsable, "Randomness in Yd");
  this->registerParam("epsilon_c",         epsilon_c,    0., _pat_parsable, "Critical strain");
  this->registerParam("Yc limit",           yc_limit, false, _pat_internal, "As the material a critical Y");
  this->registerParam("damage_in_y",   damage_in_y,   false, _pat_parsable, "Use threshold (1-D)Y");

  this->initInternalVector(this->Yd_rand, 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();

  updateInternalParameters();

  this->resizeInternalVector(this->Yd_rand);

  const Mesh & mesh = this->model->getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);

  for(; it != last_type; ++it) {
    UInt nb_element  = this->element_filter(*it).getSize();
    UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(*it);

    Vector <Real> & Yd_rand_vec = Yd_rand(*it);
    for(UInt e = 0; e < nb_element; ++e) {
      Real rand_part = (2 * drand48()-1) * Yd_randomness * Yd;

      for(UInt q = 0; q < nb_quad; ++q)
 	Yd_rand_vec(nb_quad*e+q,0) = Yd + rand_part;
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::updateInternalParameters() {
  MaterialDamage<spatial_dimension>::updateInternalParameters();
  Yc = .5 * epsilon_c * this->E * epsilon_c;
  if(std::abs(epsilon_c) > std::numeric_limits<Real>::epsilon()) yc_limit = true;
  else yc_limit = false;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::computeStress(ElementType el_type,
						      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Yd_q = Yd_rand(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Real Y = 0;

  computeStressOnQuad(grad_u, sigma, *dam, Y, *Yd_q);

  ++dam;
  ++Yd_q;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialMarigo);


__END_AKANTU__
