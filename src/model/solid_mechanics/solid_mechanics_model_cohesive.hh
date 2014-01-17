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
#include "cohesive_element_inserter.hh"
#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
#  include "facet_synchronizer.hh"
#  include "facet_stress_synchronizer.hh"
#endif
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
struct SolidMechanicsModelCohesiveOptions : public SolidMechanicsModelOptions {
  SolidMechanicsModelCohesiveOptions(AnalysisMethod analysis_method = _explicit_lumped_mass,
				     bool extrinsic = false,
				     bool no_init_materials = false,
				     bool init_facet_filter = true,
				     bool stress_interpolation = true) :
    SolidMechanicsModelOptions(analysis_method, no_init_materials),
    extrinsic(extrinsic), init_facet_filter(init_facet_filter),
    stress_interpolation(stress_interpolation) {}
  bool extrinsic;
  bool init_facet_filter;
  bool stress_interpolation;
};

extern const SolidMechanicsModelCohesiveOptions default_solid_mechanics_model_cohesive_options;

/* -------------------------------------------------------------------------- */
/* Solid Mechanics Model for Cohesive elements                                */
/* -------------------------------------------------------------------------- */

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
		const ModelOptions & options = default_solid_mechanics_model_cohesive_options);

  /// initialize the model
  void initModel();

  /// initialize cohesive material
  void initMaterials();

  /// build fragments and compute their data (mass, velocity..)
  void computeFragmentsData();

  /// init facet filters for cohesive materials
  void initFacetFilter();

  /// limit the cohesive element insertion to a given area
  void enableFacetsCheckOnArea(const Array<Real> & limits);

  /// update automatic insertion after a change in the element inserter
  void updateAutomaticInsertion();

private:

  /// initialize completely the model for extrinsic elements
  void initAutomaticInsertion();

  /// initialize stress interpolation
  void initStressInterpolation();

  /// build fragments
  void buildFragments();

  /// compute fragments' mass and velocity
  void computeFragmentsMV();

  /// compute facets' normals
  void computeNormals();

  /// resize facet stress
  void resizeFacetStress();

  /// init facets_check array
  void initFacetsCheck();

  /// fill stress_on_facet
  void fillStressOnFacet();

  /// compute average stress on elements
  void averageStressOnFacets(UInt material_index);

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */

protected:

  virtual void onNodesAdded  (const Array<UInt> & nodes_list,
			      const NewNodesEvent & event);
  virtual void onElementsAdded  (const Array<Element> & nodes_list,
				 const NewElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get facet mesh
  AKANTU_GET_MACRO(MeshFacets, mesh_facets, const Mesh &);

  /// get stress on facets vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(StressOnFacets, facet_stress, Real);

  /// get number of fragments
  AKANTU_GET_MACRO(NbFragment, nb_fragment, UInt);

  /// get mass for each fragment
  AKANTU_GET_MACRO(FragmentsMass, fragment_mass, const Array<Real> &);

  /// get average velocity for each fragment
  AKANTU_GET_MACRO(FragmentsVelocity, fragment_velocity, const Array<Real> &);

  /// get center of mass coordinates for each fragment
  AKANTU_GET_MACRO(FragmentsCenter, fragment_center, const Array<Real> &);

  /// get facet material
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(FacetMaterial, facet_material, UInt);

  /// get facet material
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(FacetMaterial, facet_material, UInt);

  /// get facet material
  AKANTU_GET_MACRO(FacetMaterial, facet_material, const ByElementTypeArray<UInt> &);

  /// THIS HAS TO BE CHANGED
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Tangents, tangents, Real);

  /// get element inserter
  AKANTU_GET_MACRO_NOT_CONST(ElementInserter, inserter, CohesiveElementInserter &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// mesh containing facets and their data structures
  Mesh & mesh_facets;

  /// @todo store tangents when normals are computed:
  ByElementTypeReal tangents;

  /// list of stresses on facet quadrature points for every element
  ByElementTypeReal stress_on_facet;

  /// stress on facets on the two sides by quadrature point
  ByElementTypeReal facet_stress;

  /// number of fragments
  UInt nb_fragment;

  /// mass for each fragment
  Array<Real> fragment_mass;

  /// average velocity for each fragment
  Array<Real> fragment_velocity;

  /// center of mass coordinates for each element
  Array<Real> fragment_center;

  /// flag to know if facets have been generated
  bool facet_generated;

  /// material to use if a cohesive element is created on a facet
  ByElementTypeUInt facet_material;

  /// stress interpolation flag
  bool stress_interpolation;

  bool is_extrinsic;

  /// cohesive element inserter
  CohesiveElementInserter inserter;

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

/* -------------------------------------------------------------------------- */
class DefaultMaterialCohesiveSelector : public DefaultMaterialSelector {
public:
  DefaultMaterialCohesiveSelector(const SolidMechanicsModelCohesive & model) :
    DefaultMaterialSelector(model.getElementIndexByMaterial()),
    facet_material(model.getFacetMaterial()),
    mesh(model.getMesh()),
    mesh_facets(model.getMeshFacets()) { }

  inline virtual UInt operator()(const Element & element) {
    if(Mesh::getKind(element.type) == _ek_cohesive) {
      try {
	const Array<Element> & cohesive_el_to_facet
	  = mesh_facets.getSubelementToElement(element.type, element.ghost_type);
	bool third_dimension = (mesh.getSpatialDimension() == 3);
	const Element & facet = cohesive_el_to_facet(element.element, third_dimension);
	return facet_material(facet.type, facet.ghost_type)(facet.element);
      } catch(...) {
	return MaterialSelector::operator()(element);
      }
    } else if (Mesh::getSpatialDimension(element.type) == mesh.getSpatialDimension() - 1) {
      return facet_material(element.type, element.ghost_type)(element.element);
    } else {
      return DefaultMaterialSelector::operator()(element);
    }
  }

private:
  const ByElementTypeUInt & facet_material;
  const Mesh & mesh;
  const Mesh & mesh_facets;
};


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SolidMechanicsModelCohesive & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__


#endif /* __AKANTU_SOLID_MECHANICS_MODEL_COHESIVE_HH__ */
