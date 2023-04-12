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
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PLANE_STRESS_TOOLBOX_HH_
#define AKANTU_PLANE_STRESS_TOOLBOX_HH_

namespace akantu {

/**
 * Empty class in dimensions different from 2
 * This class is only specialized for 2D in the tmpl file
 */
template <Int dim, class ParentMaterial = Material>
class PlaneStressToolbox : public ParentMaterial {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PlaneStressToolbox(SolidMechanicsModel & model, const ID & id = "",
                     const ID & fe_engine_id = "")
      : ParentMaterial(model, id, fe_engine_id) {}

protected:
  void initialize();

public:
  void computeAllCauchyStresses(GhostType ghost_type = _not_ghost) override {
    ParentMaterial::computeAllCauchyStresses(ghost_type);
  }

  virtual void computeCauchyStressPlaneStress(ElementType /*el_type*/,
                                              GhostType /*ghost_type*/) {
    AKANTU_DEBUG_IN();

    AKANTU_ERROR("The function \"computeCauchyStressPlaneStress\" can "
                 "only be used in 2D Plane stress problems, which means "
                 "that you made a mistake somewhere!! ");

    AKANTU_DEBUG_OUT();
  }

  virtual void computeThirdAxisDeformation(ElementType /*unused*/,
                                           GhostType /*unused*/) {}

  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    return zip_append(
        ParentMaterial::template getArguments<dim>(el_type, ghost_type),
        "C33"_n = broadcast(C33, (*this->stress)(el_type, ghost_type).size()));
  }

  decltype(auto) getArgumentsTangent(Array<Real> & tangent_matrix,
                                     ElementType el_type,
                                     GhostType ghost_type = _not_ghost) {
    return zip_append(
        ParentMaterial::template getArgumentsTangent<dim>(tangent_matrix,
                                                          el_type, ghost_type),
        "C33"_n = broadcast(C33, (*this->stress)(el_type, ghost_type).size()));
  }

protected:
  Real C33{1.};
};

#define AKANTU_PLANE_STRESS_TOOL_SPEC(dim)                                     \
  template <>                                                                  \
  inline PlaneStressToolbox<dim, Material>::PlaneStressToolbox(                \
      SolidMechanicsModel & model, const ID & id, const ID & fe_engine_id)     \
      : Material(model, id, fe_engine_id) {}

AKANTU_PLANE_STRESS_TOOL_SPEC(1)
AKANTU_PLANE_STRESS_TOOL_SPEC(3)

} // namespace akantu

#include "plane_stress_toolbox_tmpl.hh"

#endif /* AKANTU_PLANE_STRESS_TOOLBOX_HH_ */
