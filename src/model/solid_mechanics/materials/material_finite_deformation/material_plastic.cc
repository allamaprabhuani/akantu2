/**
 * @file   material_plastic.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Tue Jul 27 18:15:37 2010
 *
 * @brief  Specialization of the material class for finite deformation neo-hookean material
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
#include "material_plastic.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialPlastic<spatial_dimension>::MaterialPlastic(SolidMechanicsModel & model, const ID & id) :
Material(model, id) {
    AKANTU_DEBUG_IN();

    this->registerParam("E", E, 0., _pat_parsable | _pat_modifiable, "Young's modulus");
    this->registerParam("nu", nu, 0.5, _pat_parsable | _pat_modifiable, "Poisson's ratio");
    this->registerParam("Plane_Stress", plane_stress, false, _pat_parsmod, "Is plane stress"); /// @todo Plane_Stress should not be possible to be modified after initMaterial (but before)
    this->registerParam("lambda", lambda, _pat_readable, "First Lamé coefficient");
    this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
    this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");

    finite_deformation=true;
    //    use_previous_stress=true;

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlastic<spatial_dimension>::initMaterial() {
    AKANTU_DEBUG_IN();
    Material::initMaterial();
    if (spatial_dimension == 1) nu = 0.;
    updateInternalParameters();
    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlastic<spatial_dimension>::updateInternalParameters() {
    lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
    mu = E / (2 * (1 + nu));

    if (plane_stress) lambda = nu * E / ((1 + nu)*(1 - nu));

    kpa = lambda + 2. / 3. * mu;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialPlastic<dim>::computeStress(ElementType el_type, GhostType ghost_type) {
    AKANTU_DEBUG_IN();


    //Array<UInt> & elem_filter = element_filter(el_type, ghost_type);
    //UInt nb_element = elem_filter.getSize();

    //Array<Real>::iterator< Matrix<Real> > green_it =
    //        this->delta_strain(el_type, ghost_type).begin(dim, dim);

    //Array<Real>::iterator< Matrix<Real> > S_it =
    //        this->delta_stress(el_type, ghost_type).begin(dim, dim);

    Array<Real>::iterator< Matrix<Real> > strain_it =
            this->strain(el_type, ghost_type).begin(dim, dim);

    Array<Real>::iterator< Matrix<Real> > strain_end =
      this->strain(el_type, ghost_type).end(dim, dim);

    //Array<Real>::iterator< Matrix<Real> > previous_stress_it =
    //        this->previous_stress(el_type, ghost_type).begin(dim, dim);

    Array<Real>::iterator< Matrix<Real> > piola_it =
            this->piola_kirchhoff_stress(el_type, ghost_type).begin(dim, dim);

    //Matrix<Real> F_tensor(dim, dim);
    //Matrix<Real> E_tensor(dim, dim);
    //Matrix<Real> Piola_tensor(dim, dim);
    //Matrix<Real> invF_tensor(dim, dim);
    //Matrix<Real> invFtS(dim, dim);

    for (; strain_it != strain_end; ++strain_it, ++piola_it) {
        Matrix<Real> & grad_u = *strain_it;
        //Matrix<Real> & grad_delta_u = *green_it;
        //Matrix<Real> & delta_sigma = *S_it;
        Matrix<Real> & piola_tensor = *piola_it;
        //Matrix<Real> & cauchy_sigma = *cauchy_stress_it;

        //gradUToF<dim > (grad_u, F_tensor);
        /*switch (dim){
            case 3:
                Math::inv3(F_tensor.storage(), invF_tensor.storage());
                break;
            case 2:
                Math::inv2(F_tensor.storage(), invF_tensor.storage());
                break;

                }*/

        computeStressOnQuad(grad_u, piola_tensor);
    }

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlastic<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
        Array<Real> & tangent_matrix,
        __attribute__((unused)) GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
    computeTangentModuliOnQuad(tangent, grad_u);
    MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialPlastic<spatial_dimension>::getPushWaveSpeed() const {
    return sqrt((lambda + 2 * mu) / rho);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialPlastic<spatial_dimension>::getShearWaveSpeed() const {
    return sqrt(mu / rho);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialPlastic);

__END_AKANTU__
