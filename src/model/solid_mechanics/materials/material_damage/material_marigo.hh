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
// #include "material.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_MARIGO_HH_
#define AKANTU_MATERIAL_MARIGO_HH_

namespace akantu {

/**
 * Material marigo
 *
 * parameters in the material files :
 *   - Yd  : (default: 50)
 *   - Sd  : (default: 5000)
 *   - Ydrandomness  : (default:0)
 */
template <Int dim> class MaterialMarigo : public MaterialDamage<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using parent = MaterialDamage<dim>;

public:
  MaterialMarigo(SolidMechanicsModel & model, const ID & id = "");
  ~MaterialMarigo() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void updateInternalParameters() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

protected:
  /// constitutive law for a given quadrature point
  template <typename Args> inline void computeStressOnQuad(Args && arguments);

  template <typename Args>
  inline void computeDamageAndStressOnQuad(Args && arguments);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  [[nodiscard]] inline Int
  getNbData(const Array<Element> & elements,
            const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  decltype(auto) getArguments(ElementType el_type, GhostType ghost_type) {
    return zip_append(parent::getArguments(el_type, ghost_type),
                      "Yd"_n = this->Yd(el_type, ghost_type),
                      "Y"_n = this->Y(el_type, ghost_type));
  }

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Y, Y, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  InternalField<Real> & Y;

  /// resistance to damage
  RandomInternalField<Real> & Yd;

  /// damage threshold
  Real Sd{5000};

  /// critical epsilon when the material is considered as broken
  Real epsilon_c{0.};

  Real Yc{0.};
  bool damage_in_y{false};
  bool yc_limit{false};
};

} // namespace akantu

#include "material_marigo_inline_impl.hh"

#endif /* AKANTU_MATERIAL_MARIGO_HH_ */
