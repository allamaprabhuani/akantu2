/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_selector.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

namespace akantu {
class SolidMechanicsModelCohesive;
}

namespace akantu {

#ifndef AKANTU_MATERIAL_SELECTOR_COHESIVE_HH_
#define AKANTU_MATERIAL_SELECTOR_COHESIVE_HH_

/* -------------------------------------------------------------------------- */
/**
 * class that assigns the first cohesive material by default to the
 * cohesive elements
 */
class DefaultMaterialCohesiveSelector : public ConstitutiveLawSelector {
public:
  DefaultMaterialCohesiveSelector(const SolidMechanicsModelCohesive & model);
  Int operator()(const Element & element) override;

  AKANTU_SET_MACRO(FallbackCohesiveValue, default_cohesive_value, Int);
  AKANTU_GET_MACRO(FallbackCohesiveValue, default_cohesive_value, Int);
  static Int
  getDefaultCohesiveMaterial(const SolidMechanicsModelCohesive & model);

protected:
  const SolidMechanicsModelCohesive & model;
  const ElementTypeMapArray<Idx> & facet_material;
  const Mesh & mesh;
  Int default_cohesive_value;
};

/* -------------------------------------------------------------------------- */
/// To be used with intrinsic elements inserted along mesh physical surfaces
class MeshDataMaterialCohesiveSelector
    : public DefaultMaterialCohesiveSelector {
public:
  MeshDataMaterialCohesiveSelector(const SolidMechanicsModelCohesive & model);
  Int operator()(const Element & element) override;

private:
  const Mesh & mesh_facets;
  const ElementTypeMapArray<std::string> & material_index;
  bool third_dimension{false};
};

/// bulk1, bulk2 -> cohesive
using MaterialCohesiveRules = std::map<std::pair<ID, ID>, ID>;

/* -------------------------------------------------------------------------- */
class MaterialCohesiveRulesSelector : public DefaultMaterialCohesiveSelector {
public:
  MaterialCohesiveRulesSelector(const SolidMechanicsModelCohesive & model,
                                const MaterialCohesiveRules & rules,
                                ID mesh_data_id = "physical_names");
  Int operator()(const Element & element) override;

private:
  ID mesh_data_id;
  const Mesh & mesh_facets;
  Int spatial_dimension;
  MaterialCohesiveRules rules;
};

#endif /* AKANTU_MATERIAL_SELECTOR_COHESIVE_HH_ */

} // namespace akantu
