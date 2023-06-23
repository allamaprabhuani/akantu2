/**
 * @file   material_bulk_viscosity.hh
 *
 * @author Shenghan Zhang <shenghan.zhang@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Jul 20 2016
 *
 * @brief  Bulk viscosity in the material
 *
 * @section LICENSE
 *
 * Copyright (©) 2016 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory EESD (Earthquake Engineering and Structural Dynamics)
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_BULK_VISCOSITY_HH__
#define __AKANTU_MATERIAL_BULK_VISCOSITY_HH__

namespace akantu {

/**
 * Material bulk viscosity
 *
 * parameters in the material files :
 *   - b1  : Linear part damping coefficient (default: 0.06)
 *   - b2  : Quadratic part damping coefficient (default: 1.2)
 *
 * The default are taken from abaqus documentation, more information
 * can be found in `Numerical damping of spurious oscillations: a
 * comparison between the bulk viscosity method and the explicit
 * dissipative Tchamwa-Wielgosz scheme` L. Maheo et al., Computational
 * Mechanics (2013) 51 pages:109-128 - doi:10.1007/s00466-012-0708-8
 */
template <Int dim, class Parent = MaterialElastic<dim>>
class MaterialBulkViscosity : public Parent {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialBulkViscosity(SolidMechanicsModel & model, const ID & id = "");

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  inline decltype(auto) getArguments(ElementType el_type,
                                     GhostType ghost_type = _not_ghost) {
    return zip_append(
        Parent::getArguments(el_type, ghost_type),
        "grad_v"_n = make_view<dim, dim>(grad_v(el_type, ghost_type)),
        "l_e"_n = make_view(characteristic_length(el_type, ghost_type)),
        "sound_speed"_n = make_view(sound_speed(el_type, ghost_type)),
        "sigma_viscous"_n =
            make_view<dim, dim>(sigma_viscous(el_type, ghost_type)),
        "sigma_elastic"_n =
            make_view<dim, dim>(sigma_elastic(el_type, ghost_type)));
  }

protected:
  /// Recomputes the characteristics length and sound speed per elements
  void updateInternals();

  /// constitutive law for a given quadrature point
  template <typename Args> inline void computeStressOnQuad(Args && args) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Linear part damping coefficient
  Real b1;

  /// Quadratic part damping coefficient
  Real b2;

  /// Applied always or only on negative strain rates
  bool apply_always;

  /// Characteristic element length
  InternalField<Real> characteristic_length;

  /// sound speed
  InternalField<Real> sound_speed;

  /// Gradient of velocity
  InternalField<Real> grad_v;

  /// Elastic stress
  InternalField<Real> sigma_elastic;

  /// Viscous stress
  InternalField<Real> sigma_viscous;
};

} // namespace akantu

#include "material_bulk_viscosity_tmpl.hh"

#endif /* __AKANTU_MATERIAL_BULK_VISCOSITY_HH__ */
