/**
 * @file   solid_mechanics_model_cohesive.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Tue May 08 13:01:18 2012
 *
 * @brief  Solid mechanics model for cohesive elements
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
#include "solid_mechanics_model.hh"
#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
#  include "facet_synchronizer.hh"
#  include "facet_stress_synchronizer.hh"
#endif
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__

__BEGIN_AKANTU__

class SolidMechanicsModelCohesive : public SolidMechanicsModel {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  class NewCohesiveNodesEvent : public NewNodesEvent {
  public:
    AKANTU_GET_MACRO_NOT_CONST(OldNodesList, old_nodes, Array<UInt> &);
    AKANTU_GET_MACRO(OldNodesList, old_nodes, const Array<UInt> &);
  protected:
    Array<UInt> old_nodes;
  };

  typedef FEMTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive> MyFEMCohesiveType;

  SolidMechanicsModelCohesive(Mesh & mesh,
			      UInt spatial_dimension = _all_dimensions,
			      const ID & id = "solid_mechanics_model_cohesive",
			      const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModelCohesive();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// assemble the residual for the explicit scheme
  virtual void updateResidual(bool need_initialize = true);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// function to perform a stress check on each facet and insert
  /// cohesive elements if needed
  void checkCohesiveStress();

  /// initialize the cohesive model
  void initFull(std::string material_file,
		AnalysisMethod method = _explicit_lumped_mass,
		bool extrinsic = false,
		bool init_facet_filter = true);

  /// initialize the model
  void initModel();

  /// initialize cohesive material
  void initCohesiveMaterials(bool extrinsic,
			     bool init_facet_filter);

  /// build fragments and compute their data (mass, velocity..)
  void computeFragmentsData();

  /// init facet filters for cohesive materials
  void initFacetFilter();

  /// limit the cohesive element insertion to a given area
  void enableFacetsCheckOnArea(const Array<Real> & limits);

private:

  /// initialize completely the model for extrinsic elements
  void initAutomaticInsertion();

  /// build fragments list
  void buildFragmentsList();

  /// compute fragments' mass and velocity
  void computeFragmentsMV();

  /// compute facets' normals
  void computeNormals();

  /// resize facet array
  template<typename T>
  void resizeFacetArray(ByElementTypeArray<T> & by_el_type_array,
			bool per_quadrature_point_data = true) const;

  /// init facets_check array
  void initFacetsCheck();

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */

protected:

  virtual void onNodesAdded  (const Array<UInt> & nodes_list,
			      const NewNodesEvent & event);
  virtual void onElementsAdded  (const Array<Element> & nodes_list,
				 const NewElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get facet mesh
  AKANTU_GET_MACRO(MeshFacets, mesh_facets, const Mesh &);

  /// get facets check vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(FacetsCheck, facets_check, bool);

  /// get stress on facets vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(StressOnFacets, facet_stress, Real);

  /// get number of fragments
  AKANTU_GET_MACRO(NbFragment, nb_fragment, UInt);

  /// get fragment_to_element vectors
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(FragmentToElement, fragment_to_element, UInt);

  /// get mass for each fragment
  AKANTU_GET_MACRO(FragmentsMass, fragment_mass, const Array<Real> &);

  /// get average velocity for each fragment
  AKANTU_GET_MACRO(FragmentsVelocity, fragment_velocity, const Array<Real> &);

  /// get center of mass coordinates for each fragment
  AKANTU_GET_MACRO(FragmentsCenter, fragment_center, const Array<Real> &);

  /// get facet material
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(FacetMaterial, facet_material, UInt);

  /// THIS HAS TO BE CHANGED
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Tangents, tangents, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// mesh containing facets and their data structures
  Mesh & mesh_facets;

  /// vector containing facets in which cohesive elements can be automatically inserted
  ByElementTypeArray<bool> facets_check;

  /// @todo store tangents when normals are computed:
  ByElementTypeReal tangents;

  /// list of stresses on facet quadrature points for every element
  ByElementTypeReal stress_on_facet;

  /// stress on facets on the two sides by quadrature point
  ByElementTypeReal facet_stress;

  /// fragment number for each element
  ByElementTypeUInt fragment_to_element;

  /// number of fragments
  UInt nb_fragment;

  /// mass for each fragment
  Array<Real> fragment_mass;

  /// average velocity for each fragment
  Array<Real> fragment_velocity;

  /// center of mass coordinates for each element
  Array<Real> fragment_center;

  /// facet in which cohesive elements are inserted
  ByElementTypeArray<bool> facet_insertion;

  /// cohesive elements connected to each facet
  ByElementTypeUInt cohesive_el_to_facet;

  /// flag to know if facets have been generated
  bool facet_generated;

  /// material to use if a cohesive element is created on a facet
  ByElementTypeUInt facet_material;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
#include "solid_mechanics_model_cohesive_parallel.hh"
#endif

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_PARALLEL_COHESIVE_ELEMENT)
#  include "solid_mechanics_model_cohesive_inline_impl.cc"
#endif

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SolidMechanicsModelCohesive & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__


#endif /* __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__ */
