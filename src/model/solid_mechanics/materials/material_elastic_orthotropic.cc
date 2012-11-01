/**
 * @file   material_elastic_orthotropic.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Tue May 08 13:01:18 2012
 *
 * @brief  Orthotropic elastic material
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
#include "material_elastic_orthotropic.hh"
#include "solid_mechanics_model.hh"
#include <algorithm>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialElasticOrthotropic<spatial_dimension>::MaterialElasticOrthotropic(SolidMechanicsModel & model,
									  const ID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("E1"  , E1  , 0., _pat_parsmod, "Young's modulus (x)");
  this->registerParam("E2"  , E2  , 0., _pat_parsmod, "Young's modulus (y)");
  this->registerParam("E3"  , E3  , 0., _pat_parsmod, "Young's modulus (z)");
  this->registerParam("nu12", nu12, 0., _pat_parsmod, "Poisson's ratio (xy)");
  this->registerParam("nu13", nu13, 0., _pat_parsmod, "Poisson's ratio (xz)");
  this->registerParam("nu23", nu23, 0., _pat_parsmod, "Poisson's ratio (yz)");
  this->registerParam("G12" , G12 , 0., _pat_parsmod, "Shear modulus (xy)");
  this->registerParam("G13" , G13 , 0., _pat_parsmod, "Shear modulus (xz)");
  this->registerParam("G23" , G23 , 0., _pat_parsmod, "Shear modulus (yz)");

  plane_stress = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialElasticOrthotropic<spatial_dimension>::~MaterialElasticOrthotropic() {
  delete S;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticOrthotropic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  UInt size = this->getTangentStiffnessVoigtSize(spatial_dimension);
  S = new types::RMatrix(size, size);

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticOrthotropic<spatial_dimension>::updateInternalParameters() {
  UInt size = this->getTangentStiffnessVoigtSize(spatial_dimension);
  Real Delta;

  Real nu21 = nu12 * E2 / E1;
  Real nu31 = nu13 * E3 / E1;
  Real nu32 = nu23 * E3 / E2;

  S->clear();

  Delta = 1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 - 2 * nu21 * nu13 * nu32;
  (*S)(0,0) = E1 * (1 - nu23 * nu32) / Delta;

  if(spatial_dimension >= 2) {
    if(plane_stress) {
      Delta = 1 - nu12 * nu21;
      (*S)(0,0) = E1 / Delta;
      (*S)(1,1) = E2 / Delta;
      (*S)(0,1) = E1 * nu21 / Delta;
      (*S)(1,0) = E2 * nu12 / Delta;
      (*S)(2,2) = G12;
    }
    else {
      (*S)(0,0) = E1 * (1 - nu23 * nu32) / Delta;
      (*S)(1,1) = E2 * (1 - nu31 * nu13) / Delta;
      (*S)(0,1) = (*S)(1,0) = E2 * (nu12 + nu32 * nu13) / Delta;
      (*S)(size - 1,size - 1) = G12;
    }
    if (spatial_dimension == 3) {
      (*S)(2,2) = E3 * (1 - nu12 * nu21) / Delta;
      (*S)(0,2) = (*S)(2,0) = E1 * (nu31 + nu21 * nu32) / Delta;
      (*S)(1,2) = (*S)(2,1) = E3 * (nu23 + nu21 * nu13) / Delta;
      (*S)(4,4) = G23;
      (*S)(5,5) = G13;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticOrthotropic<spatial_dimension>::computeStress(ElementType el_type,
								  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  computeStressOnQuad(grad_u, sigma);
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticOrthotropic<spatial_dimension>::computeTangentModuli(const ElementType & el_type,
									 Vector<Real> & tangent_matrix,
									 GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  tangent.copy(*S);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElasticOrthotropic<spatial_dimension>::getPushWaveSpeed() const {
  Real Et = (*S)(0,0);
  if(spatial_dimension >= 2)
    Et = std::max(Et, (*S)(1,1));
  if(spatial_dimension == 3)
    Et = std::max(Et, (*S)(2,2));

  return sqrt( Et / rho);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElasticOrthotropic<spatial_dimension>::getShearWaveSpeed() const {
  Real G = 0;
  if(spatial_dimension == 2)
    G = (*S)(2,2);
  if(spatial_dimension == 3) {
    G = (*S)(3,3);
    G = std::max(G, (*S)(4,4));
    G = std::max(G, (*S)(5,5));
  }
  return sqrt( G / rho);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialElasticOrthotropic);


__END_AKANTU__
