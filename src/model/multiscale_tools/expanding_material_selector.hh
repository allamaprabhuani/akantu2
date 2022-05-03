/**
 * @file   expanding_material_selector.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue Jan 16 10:26:53 2014
 * @update Thursday Mar 24  2022
 * @brief  material selector for the expanding material (e.g. ASR gel)
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "aka_random_generator.hh"
#include <aka_array.hh>
#include <aka_iterators.hh>
#include <map>
#include <mesh.hh>
#include <mesh_accessor.hh>
#include <mesh_events.hh>
#include <unordered_set>
/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* ------------------------------------------------------------------ */
/* ASR material selector */
/* ------------------------------------------------------------------ */
class ExpandingMaterialSelector : public MeshDataMaterialSelector<std::string> {
public:
  ExpandingMaterialSelector(SolidMechanicsModel & model,
                            std::string expanding_material,
                            const UInt nb_expanding_elements,
                            std::string surrounding_material = "aggregate",
                            bool expanding_element_pairs = false)
      : MeshDataMaterialSelector<std::string>("physical_names", model),
        model(model), expanding_material(expanding_material),
        nb_expanding_elements(nb_expanding_elements),
        nb_placed_expanding_elements(0),
        surrounding_material(surrounding_material),
        expanding_element_pairs(expanding_element_pairs) {}

  void initExpandingElements() {
    surrounding_material_id = model.getMaterialIndex(surrounding_material);

    Mesh & mesh = this->model.getMesh();
    UInt dim = model.getSpatialDimension();
    for (auto el_type : model.getMaterial(surrounding_material)
                            .getElementFilter()
                            .elementTypes(dim)) {

      const auto & filter =
          model.getMaterial(surrounding_material).getElementFilter()(el_type);
      if (!filter.size() == 0)
        AKANTU_EXCEPTION("Check the element type for the surrounding material");

      Element el{el_type, 0, _not_ghost};
      UInt nb_element = mesh.getNbElement(el.type, el.ghost_type);
      Array<Real> barycenter(nb_element, dim);

      for (auto && data : enumerate(make_view(barycenter, dim))) {
        el.element = std::get<0>(data);
        auto & bary = std::get<1>(data);
        mesh.getBarycenter(el, bary);
      }

      /// generate the expanding elements
      UInt seed = RandomGenerator<UInt>::seed();
      std::mt19937 random_generator(seed);
      std::uniform_int_distribution<> dis(0, nb_element - 1);

      Vector<Real> center(dim);
      std::set<int> checked_baries;
      while (nb_placed_expanding_elements != nb_expanding_elements) {
        // get a random bary center
        auto bary_id = dis(random_generator);
        if (checked_baries.find(bary_id) != checked_baries.end()) {
          continue;
        }
        checked_baries.insert(bary_id);
        el.element = bary_id;
        if (MeshDataMaterialSelector<std::string>::operator()(el) ==
            surrounding_material_id) {
          expanding_elements.push_back(el);
          nb_placed_expanding_elements += 1;
        }
      }
    }
    is_expanding_mat_initialized = true;
    std::cout << nb_placed_expanding_elements << " expanding elements placed"
              << std::endl;
  }

  void initExpandingElementPairs() {
    surrounding_material_id = model.getMaterialIndex(surrounding_material);

    Mesh & mesh = this->model.getMesh();
    auto & mesh_facets = mesh.getMeshFacets();
    UInt dim = model.getSpatialDimension();
    for (auto el_type : model.getMaterial(surrounding_material)
                            .getElementFilter()
                            .elementTypes(dim)) {

      const auto & filter =
          model.getMaterial(surrounding_material).getElementFilter()(el_type);
      if (!filter.size() == 0)
        AKANTU_EXCEPTION("Check the element type for aggregate material");

      Element el{el_type, 0, _not_ghost};
      UInt nb_element = mesh.getNbElement(el.type, el.ghost_type);
      Array<Real> barycenter(nb_element, dim);

      for (auto && data : enumerate(make_view(barycenter, dim))) {
        el.element = std::get<0>(data);
        auto & bary = std::get<1>(data);
        mesh.getBarycenter(el, bary);
      }

      /// generate the expanding elements
      UInt seed = RandomGenerator<UInt>::seed();
      std::mt19937 random_generator(seed);
      std::uniform_int_distribution<> dis(0, nb_element - 1);

      Vector<Real> center(dim);
      std::set<int> checked_baries;
      while (nb_placed_expanding_elements != nb_expanding_elements) {
        // get a random bary center
        auto bary_id = dis(random_generator);
        if (checked_baries.find(bary_id) != checked_baries.end())
          continue;
        checked_baries.insert(bary_id);
        el.element = bary_id;
        if (MeshDataMaterialSelector<std::string>::operator()(el) ==
            surrounding_material_id) {
          auto & sub_to_element =
              mesh_facets.getSubelementToElement(el.type, el.ghost_type);
          auto sub_to_el_it =
              sub_to_element.begin(sub_to_element.getNbComponent());
          const Vector<Element> & subelements_to_element =
              sub_to_el_it[el.element];
          bool successful_placement{false};
          for (auto & subelement : subelements_to_element) {
            auto && connected_elements = mesh_facets.getElementToSubelement(
                subelement.type, subelement.ghost_type)(subelement.element);
            for (auto & connected_element : connected_elements) {
              if (connected_element.element == el.element)
                continue;
              if (MeshDataMaterialSelector<std::string>::operator()(
                      connected_element) == surrounding_material_id) {
                expanding_elements.push_back(el);
                expanding_elements.push_back(connected_element);
                nb_placed_expanding_elements += 1;
                successful_placement = true;
                break;
              }
            }
            if (successful_placement)
              break;
          }
        }
      }
    }
    is_expanding_mat_initialized = true;
    std::cout << nb_placed_expanding_elements
              << " expanding element pairs placed" << std::endl;
  }

  UInt operator()(const Element & elem) {

    if (not is_expanding_mat_initialized) {
      if (this->expanding_element_pairs) {
        initExpandingElementPairs();
      } else {
        initExpandingElements();
      }
    }

    UInt temp_index = MeshDataMaterialSelector<std::string>::operator()(elem);
    if (temp_index != surrounding_material_id)
      return temp_index;
    auto iit = expanding_elements.begin();
    auto eit = expanding_elements.end();
    if (std::find(iit, eit, elem) != eit) {
      return model.getMaterialIndex(expanding_material);
    }
    return temp_index;
  }

protected:
  SolidMechanicsModel & model;
  std::string expanding_material;
  std::vector<Element> expanding_elements;
  UInt nb_expanding_elements;
  UInt nb_placed_expanding_elements;
  std::string surrounding_material{"aggregate"};
  UInt surrounding_material_id{1};
  bool is_expanding_mat_initialized{false};
  bool expanding_element_pairs{false};
};

/* ------------------------------------------------------------------ */
using MaterialCohesiveRules = std::map<std::pair<ID, ID>, ID>;

class ExpandingMaterialCohesiveRulesSelector : public MaterialSelector {
public:
  ExpandingMaterialCohesiveRulesSelector(
      SolidMechanicsModelCohesive & model, const MaterialCohesiveRules & rules,
      std::string expanding_material, const UInt nb_expanding_elements,
      std::string surrounding_material = "aggregate",
      bool expanding_element_pairs = false)
      : model(model), mesh_data_id("physical_names"), mesh(model.getMesh()),
        mesh_facets(model.getMeshFacets()), dim(model.getSpatialDimension()),
        rules(rules),
        expanding_material_selector(model, expanding_material,
                                    nb_expanding_elements, surrounding_material,
                                    expanding_element_pairs),
        default_cohesive(model) {}

  UInt operator()(const Element & element) {
    if (mesh_facets.getSpatialDimension(element.type) == (dim - 1)) {
      const std::vector<Element> & element_to_subelement =
          mesh_facets.getElementToSubelement(element.type, element.ghost_type)(
              element.element);
      // Array<bool> & facets_check = model.getFacetsCheck();

      const Element & el1 = element_to_subelement[0];
      const Element & el2 = element_to_subelement[1];

      ID id1 = model.getMaterial(expanding_material_selector(el1)).getName();
      ID id2 = id1;
      if (el2 != ElementNull) {
        id2 = model.getMaterial(expanding_material_selector(el2)).getName();
      }

      auto rit = rules.find(std::make_pair(id1, id2));
      if (rit == rules.end()) {
        rit = rules.find(std::make_pair(id2, id1));
      }

      if (rit != rules.end()) {
        return model.getMaterialIndex(rit->second);
      }
    }

    if (Mesh::getKind(element.type) == _ek_cohesive) {
      return default_cohesive(element);
    }

    return expanding_material_selector(element);
  }

private:
  SolidMechanicsModelCohesive & model;
  ID mesh_data_id;
  const Mesh & mesh;
  const Mesh & mesh_facets;
  UInt dim;
  MaterialCohesiveRules rules;

  ExpandingMaterialSelector expanding_material_selector;
  DefaultMaterialCohesiveSelector default_cohesive;
};

} // namespace akantu
