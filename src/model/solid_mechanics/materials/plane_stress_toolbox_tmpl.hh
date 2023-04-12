/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "plane_stress_toolbox.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PLANE_STRESS_TOOLBOX_TMPL_HH_
#define AKANTU_PLANE_STRESS_TOOLBOX_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class ParentMaterial>
class PlaneStressToolbox<2, ParentMaterial> : public ParentMaterial {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PlaneStressToolbox(SolidMechanicsModel & model, const ID & id = "",
                     const ID & fe_engine_id = "");

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ThirdAxisDeformation,
                                         (*third_axis_deformation), Real);

protected:
  void initialize() {
    this->registerParam("Plane_Stress", plane_stress, false, _pat_parsmod,
                        "Is plane stress");
    this->third_axis_deformation =
        this->registerInternal("third_axis_deformation", 1);
    this->third_axis_deformation->setDefaultValue(1.);
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    return zip_append(
        ParentMaterial::template getArguments<2>(el_type, ghost_type),
        "C33"_n = (*this->third_axis_deformation)(el_type, ghost_type));
  }

  decltype(auto) getArgumentsTangent(Array<Real> & tangent_matrix,
                                     ElementType el_type,
                                     GhostType ghost_type = _not_ghost) {
    return zip_append(
        ParentMaterial::template getArgumentsTangent<2>(tangent_matrix, el_type,
                                                        ghost_type),
        "C33"_n = (*this->third_axis_deformation)(el_type, ghost_type));
  }

  /* ------------------------------------------------------------------------ */
  void computeStress(ElementType el_type, GhostType ghost_type) override {
    ParentMaterial::computeStress(el_type, ghost_type);
    if (this->plane_stress) {
      computeThirdAxisDeformation(el_type, ghost_type);
    }
  }

  /* ------------------------------------------------------------------------ */
  virtual void computeThirdAxisDeformation(ElementType /*unused*/,
                                           GhostType /*unused*/) {}

  /// Computation of Cauchy stress tensor in the case of finite deformation
  void computeAllCauchyStresses(GhostType ghost_type = _not_ghost) override {
    AKANTU_DEBUG_IN();

    if (this->plane_stress) {
      AKANTU_DEBUG_ASSERT(this->finite_deformation,
                          "The Cauchy stress can only be computed if you are "
                          "working in finite deformation.");

      for (auto && type : this->elementTypes(2, ghost_type)) {
        this->computeCauchyStressPlaneStress(type, ghost_type);
      }
    } else {
      ParentMaterial::computeAllCauchyStresses(ghost_type);
    }

    AKANTU_DEBUG_OUT();
  }

  virtual void
      computeCauchyStressPlaneStress(ElementType /*el_type*/,
                                     GhostType /*ghost_type*/ = _not_ghost){};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// third axis strain measure value
  std::shared_ptr<InternalField<Real>> third_axis_deformation;

  /// Plane stress or plane strain
  bool plane_stress{false};
};

template <class ParentMaterial>
inline PlaneStressToolbox<2, ParentMaterial>::PlaneStressToolbox(
    SolidMechanicsModel & model, const ID & id, const ID & fe_engine_id)
    : ParentMaterial(model, id, fe_engine_id) {
  this->initialize();
}

} // namespace akantu

#endif /* AKANTU_PLANE_STRESS_TOOLBOX_TMPL_HH_ */
