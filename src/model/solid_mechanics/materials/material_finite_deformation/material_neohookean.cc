/**
 * @file   material_neohookean.cc
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
#include "material_neohookean.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialNeohookean<spatial_dimension>::MaterialNeohookean(SolidMechanicsModel & model, const ID & id) :
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
void MaterialNeohookean<spatial_dimension>::initMaterial() {
    AKANTU_DEBUG_IN();
    Material::initMaterial();
    if (spatial_dimension == 1) nu = 0.;
    updateInternalParameters();
    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::updateInternalParameters() {
    lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
    mu = E / (2 * (1 + nu));

    if (plane_stress) lambda = nu * E / ((1 + nu)*(1 - nu));

    kpa = lambda + 2. / 3. * mu;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialNeohookean<dim>::computeStress(ElementType el_type, GhostType ghost_type) {
    AKANTU_DEBUG_IN();


    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
    computeStressOnQuad(grad_u, sigma);
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

    AKANTU_DEBUG_OUT();
}

template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computePotentialEnergy(ElementType el_type,
								GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Material::computePotentialEnergy(el_type, ghost_type);

  if(ghost_type != _not_ghost) return;
  Array<Real>::scalar_iterator epot = this->potential_energy(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computePotentialEnergyOnQuad(grad_u, *epot);
  ++epot;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
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
Real MaterialNeohookean<spatial_dimension>::getPushWaveSpeed(__attribute__((unused)) const Element & element) const {
    return sqrt((lambda + 2 * mu) / rho);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialNeohookean<spatial_dimension>::getShearWaveSpeed(__attribute__((unused)) const Element & element) const {
    return sqrt(mu / rho);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialNeohookean);

__END_AKANTU__
