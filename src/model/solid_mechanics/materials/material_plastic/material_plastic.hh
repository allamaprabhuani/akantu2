/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_PLASTIC_HH_
#define AKANTU_MATERIAL_PLASTIC_HH_

namespace akantu {

/**
 * Parent class for the plastic constitutive laws
 * parameters in the material files :
 *   - h : Hardening parameter (default: 0)
 *   - sigmay : Yield stress
 */
template <Int dim> class MaterialPlastic : public MaterialElastic<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialPlastic(SolidMechanicsModel & model, const ID & id = "",
                  const ID & fe_engine_id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * @brief Return potential or plastic energy
   *
   * Plastic dissipated energy is integrated over time.
   */
  Real getEnergy(const std::string & type) override;

  /// Update the plastic energy for the current timestep
  void updateEnergies(ElementType el_type) override;

  /// Compute the true potential energy (w/ elastic strain)
  void computePotentialEnergy(ElementType el_type) override;

protected:
  /// compute the stress and inelastic strain for the quadrature point
  template <class Args, std::enable_if_t<named_tuple_t<Args>::has(
                            "delta_grad_u"_n)> * = nullptr>
  inline void computeStressAndInelasticStrainOnQuad(Args && args) const {
    Matrix<Real, dim, dim> delta_grad_u_elastic =
        args["delta_grad_u"_n] - args["delta_inelastic_strain"_n];

    Matrix<Real, dim, dim> sigma_elastic;
    MaterialElastic<dim>::computeStressOnQuad(tuple::make_named_tuple(
        "grad_u"_n = delta_grad_u_elastic, "sigma"_n = sigma_elastic));

    args["sigma"_n] = args["previous_sigma"_n] + sigma_elastic;

    args["inelastic_strain"_n] =
        args["previous_inelastic_strain"_n] + args["delta_inelastic_strain"_n];
  }

  template <class Args, std::enable_if_t<not named_tuple_t<Args>::has(
                            "delta_grad_u"_n)> * = nullptr>
  inline void computeStressAndInelasticStrainOnQuad(Args && args) const {
    Matrix<Real, dim, dim> delta_grad_u =
        args["grad_u"_n] - args["previous_grad_u"_n];

    computeStressAndInelasticStrainOnQuad(
        tuple::append(args, "delta_grad_u"_n = delta_grad_u));
  }

  /// Get the integrated plastic energy for the time step
  Real getPlasticEnergy();

  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    return zip_append(
        MaterialElastic<dim>::getArguments(el_type, ghost_type),
        "iso_hardening"_n =
            make_view((*this->iso_hardening)(el_type, ghost_type)),
        "previous_iso_hardening"_n =
            make_view(this->iso_hardening->previous(el_type, ghost_type)),
        "inelastic_strain"_n =
            make_view<dim, dim>((*this->inelastic_strain)(el_type, ghost_type)),
        "previous_inelastic_strain"_n = make_view<dim, dim>(
            this->inelastic_strain->previous(el_type, ghost_type)));
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Yield stresss
  Real sigma_y;

  /// hardening modulus
  Real h;

  /// isotropic hardening, r
  std::shared_ptr<InternalField<Real>> iso_hardening;

  /// inelastic strain arrays ordered by element types (inelastic deformation)
  std::shared_ptr<InternalField<Real>> inelastic_strain;

  /// Plastic energy
  std::shared_ptr<InternalField<Real>> plastic_energy;

  /// @todo : add a coefficient beta that will multiply the plastic energy
  /// increment
  /// to compute the energy converted to heat

  /// Plastic energy increment
  std::shared_ptr<InternalField<Real>> d_plastic_energy;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_MATERIAL_PLASTIC_HH_ */
