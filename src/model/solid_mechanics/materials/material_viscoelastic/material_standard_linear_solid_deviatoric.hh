/**
 * @file   material_standard_linear_solid_deviatoric.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Vladislav Yastrebov <vladislav.yastrebov@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Feb 20 2020
 *
 * @brief  Material Visco-elastic, based on Standard Solid rheological model,
 * see
 * [] J.C.  Simo, T.J.R. Hughes, "Computational  Inelasticity", Springer (1998),
 * see Sections 10.2 and 10.3
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
#include "aka_common.hh"
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH_
#define AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH_

namespace akantu {

/**
 * Material standard linear solid deviatoric
 *
 *
 * @verbatim

             E_\inf
      ------|\/\/\|------
      |                 |
   ---|                 |---
      |                 |
      ----|\/\/\|--[|----
            E_v   \eta

 @endverbatim
 *
 * keyword : sls_deviatoric
 *
 * parameters in the material files :
 *   - E   : Initial Young's modulus @f$ E = E_i + E_v @f$
 *   - eta : viscosity
 *   - Ev  : stiffness of the viscous element
 */

template <Int dim>
class MaterialStandardLinearSolidDeviatoric : public MaterialElastic<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using Parent = MaterialElastic<dim>;

public:
  MaterialStandardLinearSolidDeviatoric(SolidMechanicsModel & model,
                                        const ID & id = "");
  ~MaterialStandardLinearSolidDeviatoric() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// update the internal parameters (for modifiable parameters)
  void updateInternalParameters() override;

  /// set material to steady state
  void setToSteadyState(ElementType el_type,
                        GhostType ghost_type = _not_ghost) override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  inline decltype(auto) getArguments(ElementType el_type,
                                     GhostType ghost_type = _not_ghost) {
    return zip_append(Parent::getArguments(el_type, ghost_type),
                      tuple::get<"sigma_dev"_h>() =
                          make_view<dim, dim>(stress_dev(el_type, ghost_type)),
                      tuple::get<"history"_h>() = make_view<dim, dim>(
                          history_integral(el_type, ghost_type)));
  }

protected:
  /// update the dissipated energy, is called after the stress have been
  /// computed
  void updateDissipatedEnergy(ElementType el_type, GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// give the dissipated energy for the time step
  Real getDissipatedEnergy() const;
  Real getDissipatedEnergy(const Element & element) const;

  /// get the energy using an energy type string for the time step
  Real getEnergy(const std::string & type) override;
  Real getEnergy(const std::string & energy_id,
                 const Element & element) override;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// viscosity, viscous elastic modulus
  Real eta, Ev, E_inf;

  Vector<Real> etas;

  /// history of deviatoric stress
  InternalField<Real> stress_dev;

  /// Internal variable: history integral
  InternalField<Real> history_integral;

  /// Dissipated energy
  InternalField<Real> dissipated_energy;
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH_ */
