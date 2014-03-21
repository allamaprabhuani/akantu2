/**
 * @file   material_neohookean.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 *
 * @date   Wed Jun 15 10:45:12 2011
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
#include "material_neohookean.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialNeohookean<spatial_dimension>::MaterialNeohookean(SolidMechanicsModel & model,
							  const ID & id)  :
  Material(model, id), MaterialElastic<spatial_dimension>(model, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  computeStressOnQuad(grad_u, sigma);
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;
  Real * epot = this->potential_energy(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computePotentialEnergyOnQuad(grad_u, *epot);
  epot++;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computeTangentModuli(const ElementType & el_type,
								 Array<Real> & tangent_matrix,
								 GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  computeTangentModuliOnQuad(grad_u, tangent);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialNeohookean<spatial_dimension>::celerity(const Element & elem) {
  UInt nb_quadrature_points =
    this->model->getFEM().getNbQuadraturePoints(elem.type, elem.ghost_type);

  Array<Real>::matrix_iterator strain_it = this->strain(elem.type, elem.ghost_type).begin(spatial_dimension, spatial_dimension);
  strain_it += elem.element*nb_quadrature_points;

  Real cele = 0.;

  Matrix<Real> F(3, 3);
  Matrix<Real> C(3, 3);
  Matrix<Real> Cinv(3, 3);

  for (UInt q = 0; q < nb_quadrature_points; ++q, ++strain_it) {
    Matrix<Real> & grad_u = *strain_it;

    Material::gradUToF<spatial_dimension>(grad_u, F);
    this->rightCauchy(F, C);

    Math::inv3(C.storage(), Cinv.storage());
    Real detC = Math::det3(C.storage());
    Real defvol = 0.5 * log(detC);
    Real p = this->lambda * defvol;

    Real traceCinv = Cinv.trace();
    Real coef = this->mu - p;
    Real rhot = this->rho / sqrt(detC);
    Real tmpcele = sqrt((this->lambda + 2*coef)*traceCinv*traceCinv/rhot);

    cele = std::max(cele, tmpcele);
  }

  return cele;
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialNeohookean);

__END_AKANTU__
