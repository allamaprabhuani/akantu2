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
class SolidMechanicsModel;
class FEEngine;
} // namespace akantu

namespace akantu {

/**
 * Empty class in dimensions different from 2
 * This class is only specialized for 2D in the tmpl file
 */
template <Int dim, class ParentMaterial = Material>
class AKANTU_EXPORT PlaneStressToolbox : public ParentMaterial {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PlaneStressToolbox(SolidMechanicsModel & model, const ID & id = "")
      : ParentMaterial(model, id) {}
  PlaneStressToolbox(SolidMechanicsModel & model, Int spatial_dimension,
                     const Mesh & mesh, FEEngine & fe_engine,
                     const ID & id = "")
      : ParentMaterial(model, spatial_dimension, mesh, fe_engine, id) {}

  ~PlaneStressToolbox() override = default;

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

protected:
  bool initialize_third_axis_deformation{false};
};

#define AKANTU_PLANE_STRESS_TOOL_SPEC(dim)                                     \
  template <>                                                                  \
  inline PlaneStressToolbox<dim, Material>::PlaneStressToolbox(                \
      SolidMechanicsModel & model, const ID & id)                              \
      : Material(model, id) {}

AKANTU_PLANE_STRESS_TOOL_SPEC(1)
AKANTU_PLANE_STRESS_TOOL_SPEC(3)

} // namespace akantu

#include "plane_stress_toolbox_tmpl.hh"

#endif /* AKANTU_PLANE_STRESS_TOOLBOX_HH_ */
