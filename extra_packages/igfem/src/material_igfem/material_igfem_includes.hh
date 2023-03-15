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
#ifndef AKANTU_CMAKE_LIST_MATERIALS
#include "material_igfem.hh"
#include "material_igfem_elastic.hh"
#include "material_igfem_iterative_stiffness_reduction.hh"
#include "material_igfem_saw_tooth_damage.hh"
#endif

#define AKANTU_IGFEM_MATERIAL_LIST                                             \
  ((2, (igfem_elastic, MaterialIGFEMElastic)))(                                \
      (2, (igfem_saw_tooth_damage, MaterialIGFEMSawToothDamage)))(             \
      (2, (igfem_iterative_stiffness_reduction,                                \
           MaterialIGFEMIterativeStiffnessReduction)))
