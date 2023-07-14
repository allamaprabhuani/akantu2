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
#include "material_selector_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
DefaultMaterialCohesiveSelector::DefaultMaterialCohesiveSelector(
    const SolidMechanicsModelCohesive & model)
    : model(model), facet_material(model.getFacetMaterial()),
      mesh(model.getMesh()) {
  // backward compatibility v3: to get the former behavior back when the user
  // creates its own selector
  this->setFallback(
      std::make_shared<DefaultMaterialSelector>(model.getMaterialByElement()));

  Int cohesive_index = -1;
  for (auto && material : enumerate(model.getConstitutiveLaws())) {
    if (aka::is_of_type<MaterialCohesive>(std::get<1>(material))) {
      cohesive_index = std::get<0>(material);
      break;
    }
  }

  if (cohesive_index == -1) {
    AKANTU_EXCEPTION("No cohesive materials in the material input file");
  }

  this->setFallback(cohesive_index);
}

/* -------------------------------------------------------------------------- */
Int DefaultMaterialCohesiveSelector::operator()(const Element & element) {
  if (Mesh::getKind(element.type) == _ek_cohesive) {
    try {
      const Array<Element> & cohesive_el_to_facet =
          mesh.getMeshFacets().getSubelementToElement(element.type,
                                                      element.ghost_type);
      bool third_dimension = (mesh.getSpatialDimension() == 3);
      const Element & facet =
          cohesive_el_to_facet(element.element, Int(third_dimension));
      if (facet_material.exists(facet.type, facet.ghost_type)) {
        return facet_material(facet.type, facet.ghost_type)(facet.element);
      }
      return this->getFallbackValue();

    } catch (...) {
      return this->getFallbackValue();
    }
  } else if (Mesh::getSpatialDimension(element.type) ==
             mesh.getSpatialDimension() - 1) {
    return facet_material(element.type, element.ghost_type)(element.element);
  } else {
    return ConstitutiveLawSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
MeshDataMaterialCohesiveSelector::MeshDataMaterialCohesiveSelector(
    const SolidMechanicsModelCohesive & model)
    : DefaultMaterialCohesiveSelector(model),
      mesh_facets(model.getMeshFacets()),
      material_index(mesh_facets.getData<std::string>("physical_names")),
      third_dimension(model.getSpatialDimension() == 3) {
  // backward compatibility v3: to get the former behavior back when the user
  // creates its own selector
  this->setFallback(std::make_shared<MeshDataMaterialSelector<std::string>>(
      "physical_names", model));
}

/* -------------------------------------------------------------------------- */
Int MeshDataMaterialCohesiveSelector::operator()(const Element & element) {
  if (Mesh::getKind(element.type) == _ek_cohesive or
      Mesh::getSpatialDimension(element.type) ==
          mesh_facets.getSpatialDimension() - 1) {
    Element facet{ElementNull};
    if (Mesh::getKind(element.type) == _ek_cohesive) {
      facet =
          mesh_facets.getSubelementToElement(element.type, element.ghost_type)(
              element.element, UInt(third_dimension));
    } else {
      facet = element;
    }

    try {
      std::string material_name = this->material_index(facet);
      return this->model.getMaterialIndex(material_name);
    } catch (...) {
      return this->getFallbackValue();
    }
  }
  return ConstitutiveLawSelector::operator()(element);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
MaterialCohesiveRulesSelector::MaterialCohesiveRulesSelector(
    const SolidMechanicsModelCohesive & model,
    const MaterialCohesiveRules & rules,
    ID mesh_data_id) // what we have here is the name of model and also
                     // the name of different materials
    : DefaultMaterialCohesiveSelector(model),
      mesh_data_id(std::move(mesh_data_id)), mesh_facets(model.getMeshFacets()),
      spatial_dimension(model.getSpatialDimension()), rules(rules) {

  // cohesive fallback
  this->setFallback(std::make_shared<DefaultMaterialCohesiveSelector>(model));

  // non cohesive fallback
  this->getFallbackSelector()->setFallback(
      std::make_shared<MeshDataMaterialSelector<std::string>>(
          this->mesh_data_id, model));
}

/* -------------------------------------------------------------------------- */
Int MaterialCohesiveRulesSelector::operator()(const Element & element) {
  if (mesh_facets.getSpatialDimension(element.type) ==
      (spatial_dimension - 1)) {
    const auto & element_to_subelement =
        mesh_facets.getElementToSubelement(element);
    const auto & el1 = element_to_subelement[0];
    const auto & el2 = element_to_subelement[1];

    auto id1 = this->mesh.getData<std::string>(mesh_data_id, el1);

    auto id2 = id1;
    if (el2 != ElementNull) {
      id2 = this->mesh.getData<std::string>(mesh_data_id, el2);
    }

    auto rit = rules.find(std::make_pair(id1, id2));
    if (rit == rules.end()) {
      rit = rules.find(std::make_pair(id2, id1));
    }

    if (rit != rules.end()) {
      return model.getMaterialIndex(rit->second);
    }
  }

  return ConstitutiveLawSelector::operator()(element);
}

} // namespace akantu
