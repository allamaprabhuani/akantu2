/**
 * @file   fragment_manager.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Jan 20 14:33:49 2014
 *
 * @brief  Group manager to handle fragments
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
#include "fragment_manager.hh"
#include "material_cohesive.hh"
#include <numeric>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FragmentManager::FragmentManager(SolidMechanicsModelCohesive & model,
				 bool dump_data,
				 const ID & id,
				 const MemoryID & memory_id) :
  GroupManager(model.getMesh(), id, memory_id),
  model(model),
  mass_center(0, model.getSpatialDimension(), "mass_center"),
  mass(0, model.getSpatialDimension(), "mass"),
  velocity(0, model.getSpatialDimension(), "velocity"),
  inertia_moments(0, model.getSpatialDimension(), "inertia_moments"),
  quad_coordinates("quad_coordinates", id),
  mass_density("mass_density", id),
  dump_data(dump_data) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  /// compute quadrature points' coordinates
  mesh.initByElementTypeArray(quad_coordinates,
			      spatial_dimension,
			      spatial_dimension,
			      _not_ghost);

  model.getFEM().interpolateOnQuadraturePoints(model.getMesh().getNodes(),
					       quad_coordinates);

  /// store mass density per quadrature point
  mesh.initByElementTypeArray(mass_density,
			      1,
			      spatial_dimension,
			      _not_ghost);

  storeMassDensityPerQuadraturePoint();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
class CohesiveElementFilter : public GroupManager::ClusteringFilter {
public:
  CohesiveElementFilter(const SolidMechanicsModelCohesive & model,
			const Real max_damage = 1.) :
    model(model), is_unbroken(max_damage) {}

  bool operator() (const Element & el) const {
    if (el.kind == _ek_regular)
      return true;

    const Array<UInt> & el_id_by_mat = model.getElementIndexByMaterial(el.type,
								       el.ghost_type);

    const MaterialCohesive & mat
      = static_cast<const MaterialCohesive &>
      (model.getMaterial(el_id_by_mat(el.element, 0)));

    UInt el_index = el_id_by_mat(el.element, 1);
    UInt nb_quad_per_element
      = model.getFEM("CohesiveFEM").getNbQuadraturePoints(el.type, el.ghost_type);

    const Array<Real> & damage_array = mat.getDamage(el.type, el.ghost_type);

    AKANTU_DEBUG_ASSERT(nb_quad_per_element * el_index < damage_array.getSize(),
			"This quadrature point is out of range");

    const Real * element_damage
      = damage_array.storage() + nb_quad_per_element * el_index;

    UInt unbroken_quads = std::count_if(element_damage,
					element_damage + nb_quad_per_element,
					is_unbroken);

    if (unbroken_quads > 0)
      return true;
    return false;
  }

private:

  struct IsUnbrokenFunctor {
    IsUnbrokenFunctor(const Real & max_damage) : max_damage(max_damage) {}
    bool operator() (const Real & x) {return x < max_damage;}
    const Real max_damage;
  };

  const SolidMechanicsModelCohesive & model;
  const IsUnbrokenFunctor is_unbroken;
};

/* -------------------------------------------------------------------------- */
void FragmentManager::buildFragments() {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  DistributedSynchronizer * cohesive_synchronizer
    = const_cast<DistributedSynchronizer *>(model.getCohesiveSynchronizer());

  if (cohesive_synchronizer) {
    cohesive_synchronizer->computeBufferSize(model, _gst_smmc_damage);
    cohesive_synchronizer->asynchronousSynchronize(model, _gst_smmc_damage);
    cohesive_synchronizer->waitEndSynchronize(model, _gst_smmc_damage);
  }
#endif

  DistributedSynchronizer & synchronizer
    = const_cast<DistributedSynchronizer &>(model.getSynchronizer());

  Mesh & mesh_facets = const_cast<Mesh &>(model.getMeshFacets());

  UInt spatial_dimension = model.getSpatialDimension();
  std::string fragment_prefix("fragment");

  /// generate fragments
  global_nb_fragment = createClusters(spatial_dimension,
				      fragment_prefix,
				      CohesiveElementFilter(model),
				      &synchronizer,
				      &mesh_facets);

  nb_fragment = getNbElementGroups(spatial_dimension);
  fragment_index.resize(nb_fragment);

  UInt * fragment_index_it = fragment_index.storage();

  /// loop over fragments
  for(const_element_group_iterator it(element_group_begin());
      it != element_group_end(); ++it, ++fragment_index_it) {

    /// get fragment index
    std::string fragment_index_string
      = it->first.substr(fragment_prefix.size() + 1);
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
  ByElementTypeReal unit_field("unit_field", id);
  mesh.initByElementTypeArray(unit_field,
			      spatial_dimension,
			      spatial_dimension,
			      _not_ghost);

  ByElementTypeReal::type_iterator it = unit_field.firstType(spatial_dimension,
							     _not_ghost,
							     _ek_regular);
  ByElementTypeReal::type_iterator end = unit_field.lastType(spatial_dimension,
							     _not_ghost,
							     _ek_regular);

  for (; it != end; ++it) {
    ElementType type = *it;
    Array<Real> & field_array = unit_field(type);
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_quad_per_element = model.getFEM().getNbQuadraturePoints(type);

    field_array.resize(nb_element * nb_quad_per_element);
    field_array.set(1.);
  }

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

  UInt total_components = mass_center.getSize() * mass_center.getNbComponent();

  for (UInt i = 0; i < total_components; ++i)
    mass_center_storage[i] /= mass_storage[i];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::computeVelocity() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  /// compute velocity per quadrature point
  ByElementTypeReal velocity_field("velocity_field", id);

  mesh.initByElementTypeArray(velocity_field,
			      spatial_dimension,
			      spatial_dimension,
			      _not_ghost);

  model.getFEM().interpolateOnQuadraturePoints(model.getVelocity(),
					       velocity_field);

  /// integrate on fragments
  integrateFieldOnFragments(velocity_field, velocity);

  /// divide it by the fragments' mass
  Real * mass_storage = mass.storage();
  Real * velocity_storage = velocity.storage();

  UInt total_components = velocity.getSize() * velocity.getNbComponent();

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
  ByElementTypeReal moments_coords("moments_coords", id);

  mesh.initByElementTypeArray(moments_coords,
			      spatial_dimension * spatial_dimension,
			      spatial_dimension,
			      _not_ghost);

  /// resize the by element type
  ByElementTypeReal::type_iterator it = moments_coords.firstType(spatial_dimension,
							     _not_ghost,
							     _ek_regular);
  ByElementTypeReal::type_iterator end = moments_coords.lastType(spatial_dimension,
							     _not_ghost,
							     _ek_regular);

  for (; it != end; ++it) {
    ElementType type = *it;
    Array<Real> & field_array = moments_coords(type);
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_quad_per_element = model.getFEM().getNbQuadraturePoints(type);

    field_array.resize(nb_element * nb_quad_per_element);
  }


  /// compute coordinates
  Array<Real>::const_iterator< Vector<Real> > mass_center_it
    = mass_center.begin(spatial_dimension);

  /// loop over fragments
  for(const_element_group_iterator it(element_group_begin());
      it != element_group_end(); ++it, ++mass_center_it) {

    const ByElementTypeUInt & el_list = it->second->getElements();

    ByElementTypeUInt::type_iterator type_it = el_list.firstType(spatial_dimension,
								 _not_ghost,
								 _ek_regular);
    ByElementTypeUInt::type_iterator type_end = el_list.lastType(spatial_dimension,
								 _not_ghost,
								 _ek_regular);

    /// loop over elements of the fragment
    for (; type_it != type_end; ++type_it) {
      ElementType type = *type_it;
      UInt nb_quad_per_element = model.getFEM().getNbQuadraturePoints(type);

      Array<Real> & moments_coords_array = moments_coords(type);
      const Array<Real> & quad_coordinates_array = quad_coordinates(type);
      const Array<UInt> & el_list_array = el_list(type);

      Array<Real>::iterator< Matrix<Real> > moments_begin
	= moments_coords_array.begin(spatial_dimension, spatial_dimension);
      Array<Real>::const_iterator< Vector<Real> > quad_coordinates_begin
	= quad_coordinates_array.begin(spatial_dimension);

      Vector<Real> relative_coords(spatial_dimension);

      for (UInt el = 0; el < el_list_array.getSize(); ++el) {
	UInt global_el = el_list_array(el);

	/// loop over quadrature points
	for (UInt q = 0; q < nb_quad_per_element; ++q) {
	  UInt global_q = global_el * nb_quad_per_element + q;
	  Matrix<Real> & moments_matrix = moments_begin[global_q];
	  const Vector<Real> & quad_coord_vector = quad_coordinates_begin[global_q];

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

  Array<Real>::iterator< Matrix<Real> > integrated_moments_it
    = integrated_moments.begin(spatial_dimension, spatial_dimension);
  Array<Real>::iterator< Vector<Real> > inertia_moments_it
    = inertia_moments.begin(spatial_dimension);

  for (UInt frag = 0; frag < global_nb_fragment; ++frag,
	 ++integrated_moments_it, ++inertia_moments_it)
    integrated_moments_it->eig(*inertia_moments_it);

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
void FragmentManager::storeMassDensityPerQuadraturePoint() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension);

  for (; it != end; ++it) {
    ElementType type = *it;

    Array<Real> & mass_density_array = mass_density(type);

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_quad_per_element = model.getFEM().getNbQuadraturePoints(type);
    mass_density_array.resize(nb_element * nb_quad_per_element);

    Array<UInt> & el_index_by_mat = model.getElementIndexByMaterial(type);

    Real * mass_density_it = mass_density_array.storage();

    /// store mass_density for each element and quadrature point
    for (UInt el = 0; el < nb_element; ++el) {
      Material & mat = model.getMaterial(el_index_by_mat(el, 0));

      for (UInt q = 0; q < nb_quad_per_element; ++q, ++mass_density_it)
	*mass_density_it = mat.getRho();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FragmentManager::integrateFieldOnFragments(ByElementTypeReal & field,
						Array<Real> & output) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();
  UInt nb_component = output.getNbComponent();

  /// integration part
  output.resize(global_nb_fragment);
  output.clear();

  UInt * fragment_index_it = fragment_index.storage();
  Array<Real>::iterator< Vector<Real> > output_begin = output.begin(nb_component);

  /// loop over fragments
  for(const_element_group_iterator it(element_group_begin());
      it != element_group_end(); ++it, ++fragment_index_it) {

    const ByElementTypeUInt & el_list = it->second->getElements();

    ByElementTypeUInt::type_iterator type_it = el_list.firstType(spatial_dimension,
								 _not_ghost,
								 _ek_regular);
    ByElementTypeUInt::type_iterator type_end = el_list.lastType(spatial_dimension,
								 _not_ghost,
								 _ek_regular);

    /// loop over elements of the fragment
    for (; type_it != type_end; ++type_it) {
      ElementType type = *type_it;
      UInt nb_quad_per_element = model.getFEM().getNbQuadraturePoints(type);

      const Array<Real> & density_array = mass_density(type);
      Array<Real> & field_array = field(type);
      const Array<UInt> & elements = el_list(type);
      UInt nb_element = elements.getSize();

      /// generate array to be integrated by filtering fragment's elements
      Array<Real> integration_array(elements.getSize() * nb_quad_per_element,
				    nb_component);

      Array<Real>::iterator< Matrix<Real> > int_array_it
	= integration_array.begin_reinterpret(nb_quad_per_element,
					      nb_component, nb_element);
      Array<Real>::iterator< Matrix<Real> > int_array_end
	= integration_array.end_reinterpret(nb_quad_per_element,
					    nb_component, nb_element);
      Array<Real>::iterator< Matrix<Real> > field_array_begin
	= field_array.begin_reinterpret(nb_quad_per_element,
					nb_component,
					field_array.getSize() / nb_quad_per_element);
      Array<Real>::const_iterator< Vector<Real> > density_array_begin
	= density_array.begin_reinterpret(nb_quad_per_element,
					  density_array.getSize() / nb_quad_per_element);

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
      Array<Real> integrated_array(elements.getSize(), nb_component);
      model.getFEM().integrate(integration_array,
			       integrated_array,
			       nb_component,
			       type,
			       _not_ghost,
			       elements);

      /// sum over all elements and store the result
      output_begin[*fragment_index_it]
	+= std::accumulate(integrated_array.begin(nb_component),
			   integrated_array.end(nb_component),
			   Vector<Real>(nb_component));
    }
  }

  /// sum output over all processors
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  comm.allReduce(output.storage(), global_nb_fragment * spatial_dimension, _so_sum);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void FragmentManager::createDumpDataArray(Array<T> & data,
					  std::string name,
					  bool fragment_index_output) {
  AKANTU_DEBUG_IN();

  if (data.getSize() == 0) return;

  typedef typename Array<T>::template iterator< Vector<T> > data_iterator;
  Mesh & mesh_not_const = const_cast<Mesh &>(mesh);

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_component = data.getNbComponent();
  UInt * fragment_index_it = fragment_index.storage();

  data_iterator data_begin = data.begin(nb_component);

  /// loop over fragments
  for(const_element_group_iterator it(element_group_begin());
      it != element_group_end(); ++it, ++fragment_index_it) {

    const ElementGroup & fragment = *(it->second);

    /// loop over cluster types
    ElementGroup::type_iterator type_it = fragment.firstType(spatial_dimension);
    ElementGroup::type_iterator type_end = fragment.lastType(spatial_dimension);

    for(; type_it != type_end; ++type_it) {
      ElementType type = *type_it;

      /// init mesh data
      Array<T> * mesh_data
	= mesh_not_const.getDataPointer<T>(name, type, _not_ghost, nb_component);

      data_iterator mesh_data_begin = mesh_data->begin(nb_component);

      ElementGroup::const_element_iterator el_it = fragment.element_begin(type);
      ElementGroup::const_element_iterator el_end = fragment.element_end(type);

      /// fill mesh data
      if (fragment_index_output) {
	for (; el_it != el_end; ++el_it)
	  mesh_data_begin[*el_it](0) = *fragment_index_it;
      } else {
	for (; el_it != el_end; ++el_it)
	  mesh_data_begin[*el_it] = data_begin[*fragment_index_it];
      }
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
