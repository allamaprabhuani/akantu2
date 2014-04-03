/**
 * @file   material_viscoplasticity.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 * @date   Tue Jul 09 18:15:37 20130
 *
 * @brief  Specialization of the material class for isotropic viscoplastic (small deformation) 
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
#include "material_viscoplasticity.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
MaterialViscoPlasticity<dim>::MaterialViscoPlasticity(SolidMechanicsModel & model, const ID & id) :
  Material(model, id),
  MaterialElastic<dim>(model, id),
  iso_hardening("iso_hardening", *this) {
    AKANTU_DEBUG_IN();

    this->registerParam("h", h, 0., _pat_parsable | _pat_modifiable, "Hardening  modulus");
    this->registerParam("sigmay", sigmay, 0., _pat_parsable | _pat_modifiable, "Yield stress");
    this->registerParam("rate", rate, 0., _pat_parsable | _pat_modifiable, "Rate sensitivity component");
    this->registerParam("edot0", edot0, 0., _pat_parsable | _pat_modifiable, "Reference strain rate");
    this->registerParam("ts", ts, 0., _pat_parsable | _pat_modifiable, "Time Step");
 
    this->iso_hardening.initialize(1);

    this->inelastic_deformation = true;
    this->finite_deformation    = false;
    this->use_previous_stress   = true;

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialViscoPlasticity<dim>::computeStress(ElementType el_type, GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    Array<Real>::iterator< Matrix<Real> > sigma_it =
      this->stress(el_type, ghost_type).begin(dim,dim);

    Array<Real>::iterator< Matrix<Real> > strain_it =
      this->strain(el_type, ghost_type).begin(dim,dim);

    Array<Real>::iterator< Matrix<Real> > strain_end =
      this->strain(el_type, ghost_type).end(dim,dim);

    Array<Real>::iterator< Matrix<Real> > d_strain_it =
      this->delta_strain(el_type, ghost_type).begin(dim, dim);

    Array<Real>::iterator< Matrix<Real> > inelas_strain_it =
      this->inelas_strain(el_type, ghost_type).begin(dim,dim);

    Real * iso_hardening = this->iso_hardening(el_type, ghost_type).storage();

    //Array<Real>::iterator< Matrix<Real> > previous_stress_it =
    //   this->previous_stress(el_type, ghost_type).begin(dim, dim);

    for (; strain_it != strain_end; ++strain_it, ++d_strain_it, ++sigma_it, ++inelas_strain_it, ++iso_hardening) {
      Matrix<Real> & grad_u = *strain_it;
      Matrix<Real> & grad_delta_u = *d_strain_it;
      Matrix<Real> & sigma_tensor = *sigma_it;
      Matrix<Real> & inelas_strain_tensor = *inelas_strain_it;

      computeStressOnQuad(grad_u, grad_delta_u, sigma_tensor, inelas_strain_tensor,*iso_hardening);
    }

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialViscoPlasticity<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
        Array<Real> & tangent_matrix,
        __attribute__((unused)) GhostType ghost_type) {
    AKANTU_DEBUG_IN();

        Int dim=spatial_dimension;
	//     Array<Real>::iterator< Matrix<Real> > d_strain_it =
	//       this->delta_strain(el_type, ghost_type).begin(dim, dim);

    //  Array<Real>::iterator< Matrix<Real> > sigma_it =
    //  this->stress(el_type, ghost_type).begin(dim, dim);


    Array<Real>::iterator< Matrix<Real> > previous_sigma_it =
       this->stress.previous(el_type, ghost_type).begin(dim, dim);

    Real * iso_hardening= this->iso_hardening(el_type, ghost_type).storage();


    //MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

    //   Matrix<Real> & grad_delta_u = *d_strain_it;
    //   Matrix<Real> & sigma_tensor = *sigma_it;
    //    Matrix<Real> & previous_sigma_tensor = *previous_sigma_it;

    //    computeTangentModuliOnQuad(tangent, grad_delta_u, sigma_tensor, previous_sigma_tensor, *iso_hardening);
    //   ++green_it; ++sigma_it;
    ++previous_sigma_it;
    ++iso_hardening;
    //MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

    AKANTU_DEBUG_OUT();
}

INSTANSIATE_MATERIAL(MaterialViscoPlasticity);

__END_AKANTU__
