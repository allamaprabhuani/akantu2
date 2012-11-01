/**
 * @file   material_damage_linear.cc
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the damage material
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
#include "material_damage_linear.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialDamageLinear<spatial_dimension>::MaterialDamageLinear(SolidMechanicsModel & model,
							      const ID & id)  :
  Material(model, id),
  MaterialDamage<spatial_dimension>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sigc", Sigc, 1e5, _pat_parsable, "Sigma Critique");
  this->registerParam("Gc"  , Gc  , 2. , _pat_parsable, "Gc");

  this->initInternalVector(this->K, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamageLinear<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();
  this->resizeInternalVector(this->K);
  Epsmin = Sigc / this->E;
  Epsmax = 2 * Gc/ Sigc + Epsmin;

  const Mesh & mesh = this->model->getFEM().getMesh() ;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension);
  for(; it != end; ++it) {
    UInt nb_element  = mesh.getNbElement(*it);
    UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(*it);
    Vector <Real> & K_vec = K(*it);
    std::fill_n(K_vec.storage(),nb_element*nb_quad,Epsmin);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialDamageLinear<spatial_dimension>::computeStress(ElementType el_type,
							    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * K = this->K(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  this->computeStressOnQuad(grad_u, sigma, *dam, *K);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialDamageLinear);

__END_AKANTU__
