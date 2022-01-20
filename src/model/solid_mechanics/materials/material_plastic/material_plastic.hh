/**
 * @file   material_plastic.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Common interface for plastic materials
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
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
  MaterialPlastic(SolidMechanicsModel & model, const ID & id = "");
  MaterialPlastic(SolidMechanicsModel & model, UInt a_dim, const Mesh & mesh,
                  FEEngine & fe_engine, const ID & id = "");

protected:
  void initialize();

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
  template <
      class Args,
      std::enable_if_t<tuple::has_t<"delta_grad_u"_h, Args>::value> * = nullptr>
  inline void computeStressAndInelasticStrainOnQuad(Args && args) const {
    auto && delta_grad_u_elastic = tuple::get<"delta_grad_u"_h>(args) -
                                   tuple::get<"delta_inelastic_strain"_h>(args);

    Matrix<Real, dim, dim> sigma_elastic;
    MaterialElastic<dim>::computeStressOnQuad(
        tuple::make_named_tuple(tuple::get<"grad_u"_h>() = delta_grad_u_elastic,
                                tuple::get<"sigma"_h>() = sigma_elastic));

    tuple::get<"sigma"_h>(args) =
        tuple::get<"previous_sigma"_h>(args) + sigma_elastic;

    tuple::get<"inelastic_strain"_h>(args) +=
        tuple::get<"previous_inelastic_strain"_h>(args) +
        tuple::get<"delta_inelastic_strain"_h>(args);
  }

  template <class Args, std::enable_if_t<not tuple::has_t<
                            "delta_grad_u"_h, Args>::value> * = nullptr>
  inline void computeStressAndInelasticStrainOnQuad(Args && args) const {
    Matrix<Real, dim, dim> delta_grad_u =
        tuple::get<"grad_u"_h>(args) - tuple::get<"previous_grad_u"_h>(args);

    computeStressAndInelasticStrainOnQuad(
        tuple::append(args, tuple::get<"delta_grad_u"_h>() = delta_grad_u));
  }

  /// Get the integrated plastic energy for the time step
  Real getPlasticEnergy();

  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    return zip_append(
        MaterialElastic<dim>::getArguments(el_type, ghost_type),
        tuple::get<"iso_hardening"_h>() =
            make_view(this->iso_hardening(el_type, ghost_type)),
        tuple::get<"previous_iso_hardening"_h>() =
            make_view(this->iso_hardening.previous(el_type, ghost_type)),
        tuple::get<"inelastic_strain"_h>() =
            make_view<dim, dim>(this->inelastic_strain(el_type, ghost_type)),
        tuple::get<"previous_inelastic_strain"_h>() = make_view<dim, dim>(
            this->inelastic_strain.previous(el_type, ghost_type)));
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Yield stresss
  Real sigma_y;

  /// hardening modulus
  Real h;

  /// isotropic hardening, r
  InternalField<Real> iso_hardening;

  /// inelastic strain arrays ordered by element types (inelastic deformation)
  InternalField<Real> inelastic_strain;

  /// Plastic energy
  InternalField<Real> plastic_energy;

  /// @todo : add a coefficient beta that will multiply the plastic energy
  /// increment
  /// to compute the energy converted to heat

  /// Plastic energy increment
  InternalField<Real> d_plastic_energy;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_MATERIAL_PLASTIC_HH_ */
