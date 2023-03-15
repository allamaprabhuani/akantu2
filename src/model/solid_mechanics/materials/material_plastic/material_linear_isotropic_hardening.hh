/**
 * @file   material_linear_isotropic_hardening.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Benjamin Paccaud <benjamin.paccaud@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Specialization of the material class for isotropic finite deformation
 * linear hardening plasticity
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
#include "aka_voigthelper.hh"
#include "material_plastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_LINEAR_ISOTROPIC_HARDENING_HH_
#define AKANTU_MATERIAL_LINEAR_ISOTROPIC_HARDENING_HH_

namespace akantu {

/**
 * Material plastic with a linear evolution of the yielding stress
 */
template <Int dim>
class MaterialLinearIsotropicHardening : public MaterialPlastic<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using Parent = MaterialPlastic<dim>;

public:
  using Parent::Parent;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

protected:
  /// Infinitesimal deformations
  template <class Args,
            std::enable_if_t<not named_tuple_t<Args>::has("F"_n)> * = nullptr>
  inline void computeStressOnQuad(Args && arguments);

  /// Finite deformations
  template <class Args,
            std::enable_if_t<named_tuple_t<Args>::has("F"_n)> * = nullptr>
  inline void computeStressOnQuad(Args && arguments);

  template <class Args>
  inline void computeTangentModuliOnQuad(Args && arguments) const;
};

} // namespace akantu

#include "material_linear_isotropic_hardening_inline_impl.hh"

#endif /* AKANTU_MATERIAL_LINEAR_ISOTROPIC_HARDENING_HH_ */
