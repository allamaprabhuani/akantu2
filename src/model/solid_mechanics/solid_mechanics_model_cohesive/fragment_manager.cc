/**
 * @file   fragment_manager.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Thu Jan 23 2014
 * @date last modification: Mon Dec 14 2015
 *
 * @brief  Group manager to handle fragments
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "fragment_manager.hh"
#include "aka_iterators.hh"
#include "material_cohesive.hh"
#include "mesh_iterators.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <functional>
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
FragmentManager::FragmentManager(SolidMechanicsModelCohesive & model,
                                 bool dump_data, const ID & id,
                                 const MemoryID & memory_id)
    : GroupManager(model.getMesh(), id, memory_id), model(model),
      mass_center(0, model.getSpatialDimension(), "mass_center"),
      mass(0, model.getSpatialDimension(), "mass"),
      velocity(0, model.getSpatialDimension(), "velocity"),
      inertia_moments(0, model.getSpatialDimension(), "inertia_moments"),
      principal_directions(
          0, model.getSpatialDimension() * model.getSpatialDimension(),
          "principal_directions"),
      quad_coordinates("quad_coordinates", id),
      mass_density("mass_density", id),
      nb_elements_per_fragment(0, 1, "nb_elements_per_fragment"),
      dump_data(dump_data) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  /// compute quadrature points' coordinates
  quad_coordinates.initialize(mesh, _nb_component = spatial_dimension,
                              _spatial_dimension = spatial_dimension,
                              _ghost_type = _not_ghost);
  // mesh.initElementTypeMapArray(quad_coordinates, spatial_dimension,
  //                              spatial_dimension, _not_ghost);

  model.getFEEngine().interpolateOnIntegrationPoints(model.getMesh().getNodes(),
                                                     quad_coordinates);

  /// store mass density per quadrature point
  mass_density.initialize(mesh, _spatial_dimension = spatial_dimension,
                          _ghost_type = _not_ghost);
  // mesh.initElementTypeMapArray(mass_density, 1, spatial_dimension,
  // _not_ghost);

  storeMassDensityPerIntegrationPoint();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
class CohesiveElementFilter : public GroupManager::ClusteringFilter {
public:
  CohesiveElementFilter(const SolidMechanicsModelCohesive & model,
                        const Real max_damage = 1.)
      : model(model), is_unbroken(max_damage) {}

  bool operator()(const Element & el) const {
    if (el.kind == _ek_regular)
      return true;

    const Array<UInt> & mat_indexes =
        model.getMaterialByElement(el.type, el.ghost_type);
    const Array<UInt> & mat_loc_num =
        model.getMaterialLocalNumbering(el.type, el.ghost_type);

    const MaterialCohesive & mat = static_cast<const MaterialCohesive &>(
        model.getMaterial(mat_indexes(el.element)));

    UInt el_index = mat_loc_num(el.element);
    UInt nb_quad_per_element =
        model.getFEEngine("CohesiveFEEngine")
            .getNbIntegrationPoints(el.type, el.ghost_type);

    const Array<Real> & damage_array = mat.getDamage(el.type, el.ghost_type);

    AKANTU_DEBUG_ASSERT(nb_quad_per_element * el_index < damage_array.size(),
                        "This quadrature point is out of range");

    const Real * element_damage =
        damage_array.storage() + nb_quad_per_element * el_index;

    UInt unbroken_quads = std::count_if(
        element_damage, element_damage + nb_quad_per_element, is_unbroken);

    if (unbroken_quads > 0)
      return true;
    return false;
  }

private:
  struct IsUnbrokenFunctor {
    IsUnbrokenFunctor(const Real & max_damage) : max_damage(max_damage) {}
    bool operator()(const Real & x) { return x < max_damage; }
    const Real max_damage;
  };

  const SolidMechanicsModelCohesive & model;
  const IsUnbrokenFunctor is_unbroken;
};

/* -------------------------------------------------------------------------- */
void FragmentManager::buildFragments(Real damage_limit) {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  ElementSynchronizer * cohesive_synchronizer =
      const_cast<ElementSynchronizer *>(model.getCohesiveSynchronizer());

  if (cohesive_synchronizer) {
    cohesive_synchronizer->computeBufferSize(model, _gst_smmc_damage);
    cohesive_synchronizer->asynchronousSynchronize(model, _gst_smmc_damage);
    cohesive_synchronizer->waitEndSynchronize(model, _gst_smmc_damage);
  }
#endif

  auto & mesh_facets = const_cast<Mesh &>(mesh.getMeshFacets());

  UInt spatial_dimension = model.getSpatialDimension();
  std::string fragment_prefix("fragment");

  /// generate fragments
  global_nb_fragment =
      createClusters(spatial_dimension, mesh_facets, fragment_prefix,
                     CohesiveElementFilter(model, damage_limit));

  nb_fragment = getNbElementGroups(spatial_dimension);
  fragment_index.resize(nb_fragment);

  UInt * fragment_index_it = fragment_index.storage();

  /// loop over fragments
  for (const_element_group_iterator it(element_group_begin());
       it != element_group_end(); ++it, ++fragment_index_it) {

    /// get fragment index
    std::string fragment_index_string =
        it->first.substr(fragment_prefix.size() + 1);
    std::stringstream sstr(fragment_index_string.c_str());
    sstr >> *fragment_index_it;

    AKANTU_DEBUG_ASSERT(!sstr.fail(), "fragment_index is not an integer");
  }

  /// compute fragments' mass
  computeMass();

  if (dump_data) {
    createDumpDataArray(fragment_index, "fragments", true);
    createDumpDataArray(mass, "fragments mass");
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::computeMass() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  /// create a unit field per quadrature point, since to compute mass
  /// it's neccessary to integrate only density
  ElementTypeMapArray<Real> unit_field("unit_field", id);
  unit_field.initialize(model.getFEEngine(), _nb_component = spatial_dimension,
                        _spatial_dimension = spatial_dimension,
                        _ghost_type = _not_ghost, _default_value = 1.);
  // mesh.initElementTypeMapArray(unit_field, spatial_dimension,
  // spatial_dimension,
  //                              _not_ghost);

  // ElementTypeMapArray<Real>::type_iterator it =
  //     unit_field.firstType(spatial_dimension, _not_ghost, _ek_regular);
  // ElementTypeMapArray<Real>::type_iterator end =
  //     unit_field.lastType(spatial_dimension, _not_ghost, _ek_regular);

  // for (; it != end; ++it) {
  //   ElementType type = *it;
  //   Array<Real> & field_array = unit_field(type);
  //   UInt nb_element = mesh.getNbElement(type);
  //   UInt nb_quad_per_element =
  //   model.getFEEngine().getNbIntegrationPoints(type);

  //   field_array.resize(nb_element * nb_quad_per_element);
  //   field_array.set(1.);
  // }

  integrateFieldOnFragments(unit_field, mass);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::computeCenterOfMass() {
  AKANTU_DEBUG_IN();

  /// integrate position multiplied by density
  integrateFieldOnFragments(quad_coordinates, mass_center);

  /// divide it by the fragments' mass
  Real * mass_storage = mass.storage();
  Real * mass_center_storage = mass_center.storage();

  UInt total_components = mass_center.size() * mass_center.getNbComponent();

  for (UInt i = 0; i < total_components; ++i)
    mass_center_storage[i] /= mass_storage[i];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::computeVelocity() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  /// compute velocity per quadrature point
  ElementTypeMapArray<Real> velocity_field("velocity_field", id);
  velocity_field.initialize(
      model.getFEEngine(), _nb_component = spatial_dimension,
      _spatial_dimension = spatial_dimension, _ghost_type = _not_ghost);

  // mesh.initElementTypeMapArray(velocity_field, spatial_dimension,
  //                              spatial_dimension, _not_ghost);

  model.getFEEngine().interpolateOnIntegrationPoints(model.getVelocity(),
                                                     velocity_field);

  /// integrate on fragments
  integrateFieldOnFragments(velocity_field, velocity);

  /// divide it by the fragments' mass
  Real * mass_storage = mass.storage();
  Real * velocity_storage = velocity.storage();

  UInt total_components = velocity.size() * velocity.getNbComponent();

  for (UInt i = 0; i < total_components; ++i)
    velocity_storage[i] /= mass_storage[i];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Given the distance @f$ \mathbf{r} @f$ between a quadrature point
 * and its center of mass, the moment of inertia is computed as \f[
 * I_\mathrm{CM} = \mathrm{tr}(\mathbf{r}\mathbf{r}^\mathrm{T})
 * \mathbf{I} - \mathbf{r}\mathbf{r}^\mathrm{T} \f] for more
 * information check Wikipedia
 * (http://en.wikipedia.org/wiki/Moment_of_inertia#Identities_for_a_skew-symmetric_matrix)
 *
 */

void FragmentManager::computeInertiaMoments() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  computeCenterOfMass();

  /// compute local coordinates products with respect to the center of match
  ElementTypeMapArray<Real> moments_coords("moments_coords", id);
  moments_coords.initialize(model.getFEEngine(),
                            _nb_component =
                                spatial_dimension * spatial_dimension,
                            _spatial_dimension = spatial_dimension,
                            _ghost_type = _not_ghost, _default_value = 1.);

  // mesh.initElementTypeMapArray(moments_coords,
  //                              spatial_dimension * spatial_dimension,
  //                              spatial_dimension, _not_ghost);

  // /// resize the by element type
  // ElementTypeMapArray<Real>::type_iterator it =
  //     moments_coords.firstType(spatial_dimension, _not_ghost, _ek_regular);
  // ElementTypeMapArray<Real>::type_iterator end =
  //     moments_coords.lastType(spatial_dimension, _not_ghost, _ek_regular);

  // for (; it != end; ++it) {
  //   ElementType type = *it;
  //   Array<Real> & field_array = moments_coords(type);
  //   UInt nb_element = mesh.getNbElement(type);
  //   UInt nb_quad_per_element =
  //   model.getFEEngine().getNbIntegrationPoints(type);

  //   field_array.resize(nb_element * nb_quad_per_element);
  // }

  /// compute coordinates
  Array<Real>::const_vector_iterator mass_center_it =
      mass_center.begin(spatial_dimension);

  /// loop over fragments
  for (const_element_group_iterator it(element_group_begin());
       it != element_group_end(); ++it, ++mass_center_it) {

    const ElementTypeMapArray<UInt> & el_list = it->second->getElements();
    /// loop over elements of the fragment
    for (auto type :
         el_list.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {
      UInt nb_quad_per_element =
          model.getFEEngine().getNbIntegrationPoints(type);

      Array<Real> & moments_coords_array = moments_coords(type);
      const Array<Real> & quad_coordinates_array = quad_coordinates(type);
      const Array<UInt> & el_list_array = el_list(type);

      Array<Real>::matrix_iterator moments_begin =
          moments_coords_array.begin(spatial_dimension, spatial_dimension);
      Array<Real>::const_vector_iterator quad_coordinates_begin =
          quad_coordinates_array.begin(spatial_dimension);

      Vector<Real> relative_coords(spatial_dimension);

      for (UInt el = 0; el < el_list_array.size(); ++el) {
        UInt global_el = el_list_array(el);

        /// loop over quadrature points
        for (UInt q = 0; q < nb_quad_per_element; ++q) {
          UInt global_q = global_el * nb_quad_per_element + q;
          Matrix<Real> moments_matrix = moments_begin[global_q];
          const Vector<Real> & quad_coord_vector =
              quad_coordinates_begin[global_q];

          /// to understand this read the documentation written just
          /// before this function
          relative_coords = quad_coord_vector;
          relative_coords -= *mass_center_it;

          moments_matrix.outerProduct(relative_coords, relative_coords);
          Real trace = moments_matrix.trace();
          moments_matrix *= -1.;
          moments_matrix += Matrix<Real>::eye(spatial_dimension, trace);
        }
      }
    }
  }

  /// integrate moments
  Array<Real> integrated_moments(global_nb_fragment,
                                 spatial_dimension * spatial_dimension);

  integrateFieldOnFragments(moments_coords, integrated_moments);

  /// compute and store principal moments
  inertia_moments.resize(global_nb_fragment);
  principal_directions.resize(global_nb_fragment);

  Array<Real>::matrix_iterator integrated_moments_it =
      integrated_moments.begin(spatial_dimension, spatial_dimension);
  Array<Real>::vector_iterator inertia_moments_it =
      inertia_moments.begin(spatial_dimension);
  Array<Real>::matrix_iterator principal_directions_it =
      principal_directions.begin(spatial_dimension, spatial_dimension);

  for (UInt frag = 0; frag < global_nb_fragment; ++frag,
            ++integrated_moments_it, ++inertia_moments_it,
            ++principal_directions_it) {
    integrated_moments_it->eig(*inertia_moments_it, *principal_directions_it);
  }

  if (dump_data)
    createDumpDataArray(inertia_moments, "moments of inertia");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::computeAllData() {
  AKANTU_DEBUG_IN();

  buildFragments();
  computeVelocity();
  computeInertiaMoments();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::storeMassDensityPerIntegrationPoint() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension);

  for (; it != end; ++it) {
    ElementType type = *it;

    Array<Real> & mass_density_array = mass_density(type);

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_quad_per_element = model.getFEEngine().getNbIntegrationPoints(type);
    mass_density_array.resize(nb_element * nb_quad_per_element);

    const Array<UInt> & mat_indexes = model.getMaterialByElement(type);

    Real * mass_density_it = mass_density_array.storage();

    /// store mass_density for each element and quadrature point
    for (UInt el = 0; el < nb_element; ++el) {
      Material & mat = model.getMaterial(mat_indexes(el));

      for (UInt q = 0; q < nb_quad_per_element; ++q, ++mass_density_it)
        *mass_density_it = mat.getRho();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::integrateFieldOnFragments(
    ElementTypeMapArray<Real> & field, Array<Real> & output) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();
  UInt nb_component = output.getNbComponent();

  /// integration part
  output.resize(global_nb_fragment);
  output.clear();

  UInt * fragment_index_it = fragment_index.storage();
  Array<Real>::vector_iterator output_begin = output.begin(nb_component);

  /// loop over fragments
  for (const_element_group_iterator it(element_group_begin());
       it != element_group_end(); ++it, ++fragment_index_it) {

    const ElementTypeMapArray<UInt> & el_list = it->second->getElements();

    /// loop over elements of the fragment
    for (auto type :
         el_list.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {

      UInt nb_quad_per_element =
          model.getFEEngine().getNbIntegrationPoints(type);

      const Array<Real> & density_array = mass_density(type);
      Array<Real> & field_array = field(type);
      const Array<UInt> & elements = el_list(type);
      UInt nb_element = elements.size();

      /// generate array to be integrated by filtering fragment's elements
      Array<Real> integration_array(elements.size() * nb_quad_per_element,
                                    nb_component);

      Array<Real>::matrix_iterator int_array_it =
          integration_array.begin_reinterpret(nb_quad_per_element, nb_component,
                                              nb_element);
      Array<Real>::matrix_iterator int_array_end =
          integration_array.end_reinterpret(nb_quad_per_element, nb_component,
                                            nb_element);
      Array<Real>::matrix_iterator field_array_begin =
          field_array.begin_reinterpret(nb_quad_per_element, nb_component,
                                        field_array.size() /
                                            nb_quad_per_element);
      Array<Real>::const_vector_iterator density_array_begin =
          density_array.begin_reinterpret(nb_quad_per_element,
                                          density_array.size() /
                                              nb_quad_per_element);

      for (UInt el = 0; int_array_it != int_array_end; ++int_array_it, ++el) {
        UInt global_el = elements(el);
        *int_array_it = field_array_begin[global_el];

        /// multiply field by density
        const Vector<Real> & density_vector = density_array_begin[global_el];

        for (UInt i = 0; i < nb_quad_per_element; ++i) {
          for (UInt j = 0; j < nb_component; ++j) {
            (*int_array_it)(i, j) *= density_vector(i);
          }
        }
      }

      /// integrate the field over the fragment
      Array<Real> integrated_array(elements.size(), nb_component);
      model.getFEEngine().integrate(integration_array, integrated_array,
                                    nb_component, type, _not_ghost, elements);

      /// sum over all elements and store the result
      Vector<Real> output_tmp(output_begin[*fragment_index_it]);
      output_tmp += std::accumulate(integrated_array.begin(nb_component),
                                    integrated_array.end(nb_component),
                                    Vector<Real>(nb_component));
    }
  }

  /// sum output over all processors
  const StaticCommunicator & comm = mesh.getCommunicator();
  comm.allReduce(output, _so_sum);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::computeNbElementsPerFragment() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();
  nb_elements_per_fragment.resize(global_nb_fragment);
  nb_elements_per_fragment.clear();

  UInt * fragment_index_it = fragment_index.storage();

  /// loop over fragments
  for (const_element_group_iterator it(element_group_begin());
       it != element_group_end(); ++it, ++fragment_index_it) {

    const ElementTypeMapArray<UInt> & el_list = it->second->getElements();

    /// loop over elements of the fragment
    for (auto type :
         el_list.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {
      UInt nb_element = el_list(type).size();

      nb_elements_per_fragment(*fragment_index_it) += nb_element;
    }
  }

  /// sum values over all processors
  const auto & comm = mesh.getCommunicator();
  comm.allReduce(nb_elements_per_fragment, _so_sum);

  if (dump_data)
    createDumpDataArray(nb_elements_per_fragment, "elements per fragment");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void FragmentManager::createDumpDataArray(Array<T> & data, std::string name,
                                          bool fragment_index_output) {
  AKANTU_DEBUG_IN();

  if (data.size() == 0)
    return;

  auto & mesh_not_const = const_cast<Mesh &>(mesh);

  auto && spatial_dimension = mesh.getSpatialDimension();
  auto && nb_component = data.getNbComponent();
  auto && data_begin = data.begin(nb_component);
  auto fragment_index_it = fragment_index.begin();

  /// loop over fragments
  for (const auto & fragment : ElementGroupsIterable(*this)) {
    const auto & fragment_idx = *fragment_index_it;

    /// loop over cluster types
    for (auto & type : fragment.elementTypes(spatial_dimension)) {
      /// init mesh data
      auto & mesh_data = mesh_not_const.getDataPointer<T>(
          name, type, _not_ghost, nb_component);

      auto mesh_data_begin = mesh_data.begin(nb_component);

      /// fill mesh data
      for (const auto & elem : fragment.getElements(type)) {
        Vector<T> md_tmp = mesh_data_begin[elem];
        if (fragment_index_output) {
          md_tmp(0) = fragment_idx;
        } else {
          md_tmp = data_begin[fragment_idx];
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
