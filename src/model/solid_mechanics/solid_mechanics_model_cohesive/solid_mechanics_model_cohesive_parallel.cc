/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communicator.hh"
#include "element_synchronizer.hh"
#include "material_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <type_traits>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
// void SolidMechanicsModelCohesive::synchronizeGhostFacetsConnectivity() {
//   AKANTU_DEBUG_IN();

//   const Communicator & comm = mesh.getCommunicator();
//   Int psize = comm.getNbProc();

//   if (psize == 1) {
//     AKANTU_DEBUG_OUT();
//     return;
//   }

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateCohesiveSynchronizers(
    NewElementsEvent & elements_event) {
  /// update synchronizers if needed

  if (not mesh.isDistributed()) {
    return;
  }

  ElementTypeMap<Int> nb_new_cohesive_elements;
  for (auto ghost_type : ghost_types) {
    for (auto cohesive_type :
         mesh.elementTypes(spatial_dimension, ghost_type, _ek_cohesive)) {
      nb_new_cohesive_elements(cohesive_type, ghost_type) = 0;
    }
  }

  for (auto & el : elements_event.getList()) {
    if (el.kind() != _ek_cohesive) {
      continue;
    }
    ++nb_new_cohesive_elements(el.type, el.ghost_type);
  }

  auto & mesh_facets = inserter->getMeshFacets();
  auto & facet_synchronizer = mesh_facets.getElementSynchronizer();
  const auto & cfacet_synchronizer = facet_synchronizer;

  // update the cohesive element synchronizer
  cohesive_synchronizer->updateSchemes(
      [&](auto && scheme, auto && proc, auto && direction) {
        auto & facet_scheme =
            cfacet_synchronizer.getCommunications().getScheme(proc, direction);

        for (auto && facet : facet_scheme) {
          const auto & cohesive_element =
              mesh_facets.getElementToSubelement(facet)[1];

          if (cohesive_element == ElementNull or
              cohesive_element.kind() != _ek_cohesive) {
            continue;
          }

          auto && cohesive_type = FEEngine::getCohesiveElementType(facet.type);
          auto old_nb_cohesive_elements =
              mesh.getNbElement(cohesive_type, facet.ghost_type);
          old_nb_cohesive_elements -=
              nb_new_cohesive_elements(cohesive_type, facet.ghost_type);

          if (cohesive_element.element >= old_nb_cohesive_elements) {
            scheme.push_back(cohesive_element);
          }
        }
      });

  if (not facet_stress_synchronizer) {
    return;
  }

  const auto & element_synchronizer = mesh.getElementSynchronizer();
  const auto & comm = mesh.getCommunicator();
  auto && my_rank = comm.whoAmI();

  // update the facet stress synchronizer
  facet_stress_synchronizer->updateSchemes([&](auto && scheme, auto && proc,
                                               auto && /*direction*/) {
    auto it_element = scheme.begin();
    for (auto && element : scheme) {
      auto && facet_check = inserter->getCheckFacets(
          element.type, element.ghost_type)(element.element); // slow access
                                                              // here

      if (facet_check) {
        auto && connected_elements = mesh_facets.getElementToSubelement(
            element.type, element.ghost_type)(element.element); // slow access
                                                                // here
        auto && rank_left = element_synchronizer.getRank(connected_elements[0]);
        auto && rank_right =
            element_synchronizer.getRank(connected_elements[1]);

        // keep element if the element is still a boundary element between two
        // processors
        if ((rank_left == Int(proc) and rank_right == my_rank) or
            (rank_left == my_rank and rank_right == Int(proc))) {
          *it_element = element;
          ++it_element;
        }
      }
    }
    scheme.resize(it_element - scheme.begin());
  });
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateFacetStressSynchronizer() {
  if (facet_stress_synchronizer == nullptr) {
    return;
  }

  const auto & rank_to_element =
      mesh.getElementSynchronizer().getElementToRank();
  const auto & facet_checks = inserter->getCheckFacets();
  const auto & mesh_facets = inserter->getMeshFacets();
  const auto & element_to_facet = mesh_facets.getElementToSubelement();
  auto rank = mesh.getCommunicator().whoAmI();

  facet_stress_synchronizer->updateSchemes(
      [&](auto & scheme, auto & proc, auto & /*direction*/) {
        Idx el = 0;
        for (auto && element : scheme) {
          if (not facet_checks(element)) {
            continue;
          }

          const auto & next_el = element_to_facet(element);
          auto rank_left = rank_to_element(next_el[0]);
          auto rank_right = rank_to_element(next_el[1]);

          if ((rank_left == rank and rank_right == proc) or
              (rank_left == proc and rank_right == rank)) {
            scheme[el] = element;
            ++el;
          }
        }
        scheme.resize(el);
      });
}

/* -------------------------------------------------------------------------- */
template <typename T>
void SolidMechanicsModelCohesive::packFacetStressDataHelper(
    const ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & elements) const {
  packUnpackFacetStressDataHelper<T, true>(data_to_pack, buffer, elements);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void SolidMechanicsModelCohesive::unpackFacetStressDataHelper(
    ElementTypeMapArray<T> & data_to_unpack, CommunicationBuffer & buffer,
    const Array<Element> & elements) const {
  packUnpackFacetStressDataHelper<T, false>(data_to_unpack, buffer, elements);
}

/* -------------------------------------------------------------------------- */
template <typename T, bool pack_helper>
void SolidMechanicsModelCohesive::packUnpackFacetStressDataHelper(
    std::conditional_t<pack_helper, const ElementTypeMapArray<T>,
                       ElementTypeMapArray<T>> & data_to_pack,
    CommunicationBuffer & buffer, const Array<Element> & elements) const {
  bool element_rank = false;
  auto & mesh_facets = inserter->getMeshFacets();

  auto & fe_engine = this->getFEEngine("FacetsFEEngine");
  for (auto && el : elements) {
    if (el.type == _not_defined) {
      AKANTU_EXCEPTION(
          "packUnpackFacetStressDataHelper called with wrong inputs");
    }

    auto ghost_type = mesh_facets.getElementToSubelement()(el)[0].ghost_type;
    if constexpr (pack_helper) {
      element_rank = ghost_type != _not_ghost;
    } else {
      element_rank = ghost_type == _not_ghost;
    }

    auto nb_quad_per_elem =
        fe_engine.getNbIntegrationPoints(el.type, el.ghost_type);

    auto && data_per_element = data_to_pack.get(
        el, spatial_dimension * spatial_dimension, 2, nb_quad_per_elem);
    for (auto && data_per_quad : data_per_element) {
      auto && data = data_per_quad(element_rank);

      if constexpr (pack_helper) {
        buffer << data;
      } else {
        buffer >> data;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
Int SolidMechanicsModelCohesive::getNbQuadsForFacetCheck(
    const Array<Element> & elements) const {
  Int nb_quads = 0;
  const auto & fe_engine = this->getFEEngine("FacetsFEEngine");
  for (const auto & el : elements) {
    auto nb_quad_per_facet =
        fe_engine.getNbIntegrationPoints(el.type, el.ghost_type);

    nb_quads += nb_quad_per_facet;
  }

  return nb_quads;
}

/* -------------------------------------------------------------------------- */
Int SolidMechanicsModelCohesive::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  Int size = 0;
  if (elements.empty()) {
    return 0;
  }

  /// regular element case
  if (elements(0).kind() == _ek_regular) {
    switch (tag) {
    // case SynchronizationTag::_smmc_facets: {
    //   size += elements.size() * sizeof(bool);
    //   break;
    // }
    case SynchronizationTag::_smmc_facets_stress: {
      auto nb_quads = getNbQuadsForFacetCheck(elements);
      size += nb_quads * spatial_dimension * spatial_dimension * sizeof(Real);
      break;
    }
    case SynchronizationTag::_constitutive_law_id: {
      for (auto && element : elements) {
        if (Mesh::getSpatialDimension(element.type) ==
            (spatial_dimension - 1)) {
          size += sizeof(Idx);
        }
      }

      size += SolidMechanicsModel::getNbData(elements, tag);
      break;
    }

    default: {
      size += SolidMechanicsModel::getNbData(elements, tag);
    }
    }
  } else if (elements(0).kind() == _ek_cohesive) {

    switch (tag) {
    case SynchronizationTag::_smm_boundary: {
      Int nb_nodes_per_element = 0;

      for (auto && el : elements) {
        nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
      }

      // force, displacement, boundary
      size += nb_nodes_per_element * spatial_dimension *
              Int(2 * sizeof(Real) + sizeof(bool));
      break;
    }
    default:
      break;
    }

    if (tag != SynchronizationTag::_smmc_facets) {
      size += CLHParent::getNbData(elements, tag);
    }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  if (elements.empty()) {
    return;
  }

  if (elements(0).kind() == _ek_regular) {
    switch (tag) {
    case SynchronizationTag::_smmc_facets_stress: {
      packFacetStressDataHelper(facet_stress, buffer, elements);
      break;
    }
    case SynchronizationTag::_constitutive_law_id: {
      for (auto && element : elements) {
        if (Mesh::getSpatialDimension(element.type) !=
            (spatial_dimension - 1)) {
          continue;
        }
        buffer << this->getConstitutiveLawByElement()(element);
      }

      SolidMechanicsModel::packData(buffer, elements, tag);
      break;
    }
    default: {
      SolidMechanicsModel::packData(buffer, elements, tag);
    }
    }

    AKANTU_DEBUG_OUT();
    return;
  }

  if (elements(0).kind() == _ek_cohesive) {
    switch (tag) {
    case SynchronizationTag::_smm_boundary: {
      packNodalDataHelper(*internal_force, buffer, elements, mesh);
      packNodalDataHelper(*velocity, buffer, elements, mesh);
      packNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
      break;
    }
    default: {
    }
    }

    if (tag != SynchronizationTag::_smmc_facets) {
      CLHParent::packData(buffer, elements, tag);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::unpackData(CommunicationBuffer & buffer,
                                             const Array<Element> & elements,
                                             const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (elements.empty()) {
    return;
  }

  if (elements(0).kind() == _ek_regular) {
    switch (tag) {
    case SynchronizationTag::_smmc_facets_stress: {
      unpackFacetStressDataHelper(facet_stress, buffer, elements);
      break;
    }
    case SynchronizationTag::_constitutive_law_id: {
      for (auto && element : elements) {
        if (Mesh::getSpatialDimension(element.type) !=
            (spatial_dimension - 1)) {
          continue;
        }

        Int recv_mat_index;
        buffer >> recv_mat_index;
        auto & mat_index = getConstitutiveLawByElement()(element);
        if (mat_index != Int(-1)) {
          continue;
        }

        // add ghosts element to the correct material
        mat_index = recv_mat_index;
        auto & mat =
            aka::as_type<MaterialCohesive>(getConstitutiveLaw(mat_index));
        if (is_extrinsic) {
          mat.addFacet(element);
        }
        facet_material(element) = recv_mat_index;
      }
      SolidMechanicsModel::unpackData(buffer, elements, tag);
      break;
    }
    default: {
      SolidMechanicsModel::unpackData(buffer, elements, tag);
    }
    }

    AKANTU_DEBUG_OUT();
    return;
  }

  if (elements(0).kind() == _ek_cohesive) {
    switch (tag) {
    case SynchronizationTag::_smm_boundary: {
      unpackNodalDataHelper(*internal_force, buffer, elements, mesh);
      unpackNodalDataHelper(*velocity, buffer, elements, mesh);
      unpackNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
      break;
    }
    default: {
    }
    }

    if (tag != SynchronizationTag::_smmc_facets) {
      CLHParent::unpackData(buffer, elements, tag);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
