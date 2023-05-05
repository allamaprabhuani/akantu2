/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AKANTU_MATERIAL_EXTRA_INCLUDES_HH_
#define AKANTU_MATERIAL_EXTRA_INCLUDES_HH_

#ifndef AKANTU_CMAKE_LIST_MATERIALS

// visco-elastic materials
#include "material_stiffness_proportional.hh"

// damage materials
#include "material_brittle.hh"
#include "material_damage_iterative.hh"
#include "material_damage_linear.hh"
#include "material_iterative_stiffness_reduction.hh"
#include "material_orthotropic_damage_iterative.hh"
#include "material_vreepeerlings.hh"

// plasticity
#include "material_viscoplastic.hh"

// multi-scale simulations
#include "material_FE2.hh"

#endif

#if defined(AKANTU_DAMAGE_NON_LOCAL)
#ifndef AKANTU_CMAKE_LIST_MATERIALS
#include "material_brittle_non_local.hh"
#include "material_damage_iterative_non_local.hh"
#include "material_orthotropic_damage_iterative_non_local.hh"
#include "material_orthotropic_damage_non_local.hh"
#include "material_vreepeerlings_non_local.hh"
#endif

#define AKANTU_DAMAGE_NON_LOCAL_MATERIAL_EXTRA_LIST                            \
  ((2, (brittle_non_local, MaterialBrittleNonLocal)))(                         \
      (2, (damage_iterative_non_local, MaterialDamageIterativeNonLocal)))(     \
      (2, (damage_orthotropic_iterative_non_local,                             \
           MaterialOrthotropicDamageIterativeNonLocal)))
#else
#define AKANTU_DAMAGE_NON_LOCAL_MATERIAL_EXTRA_LIST
#endif

#define AKANTU_EXTRA_MATERIAL_LIST                                             \
  ((2, (damage_linear, MaterialDamageLinear)))(                                \
      (2, (brittle, MaterialBrittle)))((2, (material_FE2, MaterialFE2)))(      \
      (2, (damage_iterative, MaterialDamageIterative)))(                       \
      (2,                                                                      \
       (iterative_stiffness_reduction, MaterialIterativeStiffnessReduction)))( \
      (2, (vreepeerlings, MaterialVreePeerlings)))(                            \
      (2, (ve_stiffness_prop, MaterialStiffnessProportional)))(                \
      (2, (visco_plastic, MaterialViscoPlastic)))(                             \
      (2, (orthotropic_damage_iterative, MaterialOrthotropicDamageIterative)))

#endif /* AKANTU_MATERIAL_EXTRA_INCLUDES_HH_ */
