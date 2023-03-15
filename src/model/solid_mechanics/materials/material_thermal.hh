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
#include "aka_common.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_THERMAL_HH_
#define AKANTU_MATERIAL_THERMAL_HH_

namespace akantu {
template <Int dim> class MaterialThermal : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialThermal(SolidMechanicsModel & model, const ID & id = "");
  MaterialThermal(SolidMechanicsModel & model, Int spatial_dimension,
                  const Mesh & mesh, FEEngine & fe_engine, const ID & id = "");

  ~MaterialThermal() override = default;

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type) override;

  /// local computation of thermal stress
  template <class Args> inline void computeStressOnQuad(Args && args);

  /* ------------------------------------------------------------------------ */
  decltype(auto) getArguments(ElementType el_type, GhostType ghost_type) {
    return zip_append(
        Material::getArguments<dim>(el_type, ghost_type),
        "delta_T"_n = make_view(this->delta_T(el_type, ghost_type)),
        "sigma_th"_n = make_view(this->sigma_th(el_type, ghost_type)),
        "previous_sigma_th"_n =
            make_view(this->sigma_th.previous(el_type, ghost_type)));
  }

  decltype(auto) getArgumentsTangent(Array<Real> & tangent_matrices,
                                     ElementType el_type,
                                     GhostType ghost_type) {
    return Material::getArgumentsTangent<dim>(tangent_matrices, el_type,
                                              ghost_type);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Young modulus
  Real E;

  /// Poisson ratio
  Real nu;

  /// Thermal expansion coefficient
  /// TODO : implement alpha as a matrix
  Real alpha;

  /// Temperature field
  InternalField<Real> delta_T;

  /// Current thermal stress
  InternalField<Real> sigma_th;
};

/* ------------------------------------------------------------------------ */
/* Inline impl                                                              */
/* ------------------------------------------------------------------------ */
template <Int dim>
template <class Args>
inline void MaterialThermal<dim>::computeStressOnQuad(Args && args) {
  auto && sigma = args["sigma_th"_n];
  auto && deltaT = args["delta_T"_n];
  sigma = -this->E / (1. - 2. * this->nu) * this->alpha * deltaT;
}

template <>
template <class Args>
inline void MaterialThermal<1>::computeStressOnQuad(Args && args) {
  auto && sigma = args["sigma_th"_n];
  auto && deltaT = args["delta_T"_n];
  sigma = -this->E * this->alpha * deltaT;
}

} // namespace akantu

#endif /* AKANTU_MATERIAL_THERMAL_HH_ */
