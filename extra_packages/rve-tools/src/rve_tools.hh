/**
 * @file   rve_tools.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue Jan 16 10:26:53 2014
 * @update Tue Feb 8  2022
 * @brief  tools for the analysis of multiscale problems and single RVEs
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
#include <boundary_condition_functor.hh>
#include <map>
#include <mesh.hh>
#include <mesh_accessor.hh>
#include <mesh_events.hh>
#include <unordered_set>
/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear_sequential.hh"
#include "solid_mechanics_model.hh"
#include "solid_mechanics_model_cohesive.hh"
#ifdef AKANTU_FLUID_DIFFUSION
#include "fluid_diffusion_model.hh"
#endif
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_RVE_TOOLS_HH__
#define __AKANTU_RVE_TOOLS_HH__

namespace akantu {
class NodeGroup;
class SolidMechanicsModel;
} // namespace akantu

namespace akantu {

class RVETools : public MeshEventHandler {
public:
  RVETools(SolidMechanicsModel & model);

  virtual ~RVETools() = default;

  /// This function is used in case of uni-axial boundary conditions:
  /// It detects the nodes, on which a traction needs to be applied and stores
  /// them in a node group
  void fillNodeGroup(NodeGroup & node_group, SpatialDirection dir,
                     Real offset = 0);

  /// compute volumes of each phase
  void computePhaseVolumes();

  /// compute volume of the model passed to RVETools (RVE if FE2_mat)
  void computeModelVolume();

  /// Apply free expansion boundary conditions
  void applyFreeExpansionBC(Real offset = 0);

  /// Apply loaded boundary conditions
  void applyLoadedBC(const Vector<Real> & traction, const ID & element_group,
                     bool multi_axial);
  /// Apply traction on the upper limit in the given direction
  void applyExternalTraction(Real traction, SpatialDirection direction);

  /// This function computes the average strain along one side
  Real computeAverageStrain(SpatialDirection direction, Real offset = 0);

  /// Computes a compoment of average volumetric strain tensor,
  /// i.e. eps_macro_xx or eps_macro_yy or eps_macro_zz
  Real computeVolumetricExpansion(SpatialDirection direction);

  /// This function computes the amount of damage in one material phase
  Real computeDamagedVolume(const ID & mat_name);

  /// This function is used to compute the average stiffness by performing a
  /// virtual loading test
  Real performLoadingTest(SpatialDirection direction, bool tension);

  /// perform tension tests and integrate the internal force on the upper
  /// surface
  void computeStiffnessReduction(std::ofstream & file_output, Real time,
                                 bool tension = true);

  /// just the average properties, NO tension test
  void computeAverageProperties(std::ofstream & file_output);

  /// Same as the last one but writes a csv with first column having time in
  /// days
  void computeAverageProperties(std::ofstream & file_output, Real time);

  /// computes and writes only displacements and strains
  void computeAveragePropertiesCohesiveModel(std::ofstream & file_output,
                                             Real time, Real offset = 0);

  /// Averages strains & displacements + damage ratios in FE2 material
  void computeAveragePropertiesFe2Material(std::ofstream & file_output,
                                           Real time);

  /// Save the state of the simulation for subsequent restart
  void saveState(UInt current_step);

  /// Load to previous state of the simulation for restart
  bool loadState(UInt & current_step);

  /// compute the volume occupied by a given material
  Real computePhaseVolume(const ID & mat_name);

  /// apply eigenstrain in cracks
  void applyEigenGradUinCracks(const Matrix<Real> prescribed_eigen_grad_u,
                               const ID & material_name);

  /// apply eigenstrain in the cracks, that advanced in the last step
  void applyEigenGradUinCracks(const Matrix<Real> prescribed_eigen_grad_u,
                               const ElementTypeMapUInt & critical_elements,
                               bool to_add = false);

  Real computeSmallestElementSize();

  // /// apply homogeneous temperature on Solid Mechanics Model
  // void applyTemperatureFieldToSolidmechanicsModel(const Real & temperature);

  /// compute chemical expandion by a sigmoidal rule (Larive, 1998)
  void computeChemExpansionLarive(const Real & delta_time_day, const Real & T,
                                  Real & expansion, const Real & eps_inf,
                                  const Real & time_ch_ref,
                                  const Real & time_lat_ref, const Real & U_C,
                                  const Real & U_L, const Real & T_ref);

  /// compute chemical expansion based on the Arrhenius equation
  void computeChemExpansionArrhenius(const Real & delta_time_day,
                                     const Real & T, Real & expansion,
                                     const Real & k, const Real & Ea);

  /// compute increase in eigen strain within 1 timestep
  Real computeDeltaEigenStrainThermal(const Real delta_time_day, const Real k,
                                      const Real activ_energy, const Real R,
                                      const Real T,
                                      Real & amount_reactive_particles,
                                      const Real saturation_const);

  /// compute linear increase in eigen strain
  Real computeDeltaEigenStrainLinear(const Real delta_time, const Real k);

  /// insert multiple blocks of cohesive elements
  void insertCohesivesRandomly(const UInt & nb_coh_elem,
                               std::string matrix_mat_name, Real gap_ratio);
  /// insert multiple blocks of cohesive elements
  void insertCohesivesRandomly3D(const UInt & nb_coh_elem,
                                 std::string matrix_mat_name, Real gap_ratio);

  /// crack insertion for 1st order 3D elements
  void insertCohesiveLoops3D(const UInt & nb_insertions,
                             std::string facet_mat_name, Real gap_ratio);

  /// insert pre cracks and output number of successful insertions
  template <UInt dim>
  UInt insertPreCracks(const UInt & nb_insertions, std::string facet_mat_name);

  /// insert block of cohesive elements based on the coord of the central
  void insertCohesivesByCoords(const Matrix<Real> & positions,
                               Real gap_ratio = 0);

  /// communicates crack numbers from not ghosts to ghosts cohesive elements
  void communicateCrackNumbers();

protected:
  /// put flags on all nodes who have a ghost counterpart
  void communicateFlagsOnNodes();

  /// on master or slave nodes with potential facet loops around
  /// put 0 eff stress for the smallest value
  void communicateEffStressesOnNodes();

  /// pick facets by passed coordinates
  void pickFacetsByCoord(const Matrix<Real> & positions);

  /// pick facets randomly within specified material
  void pickFacetsRandomly(UInt nb_insertions, std::string facet_mat_name);

  /// build closed facets loop around random point of specified material and
  /// returns number of successfully inserted loops
  UInt closedFacetsLoopAroundPoint(UInt nb_insertions, std::string mat_name);

  /// finds facet loops around a point and insert cohesives
  template <UInt dim>
  UInt insertCohesiveLoops(UInt nb_insertions, std::string mat_name);

  /// builds a facet loop around a point from starting to ending segm-s
  Array<Element> findFacetsLoopFromSegment2Segment(
      Element starting_facet, Element starting_segment, Element ending_segment,
      UInt cent_node, Real max_dot, UInt material_id, bool check_crack_facets);

  /// builds a facet loop around a point from starting to ending segment using
  /// the Dijkstra shortest path algorithm in boost library
  /// uses sum of distances to 2 incenters as weight
  Array<Element>
  findFacetsLoopByGraphByDist(const Array<Element> & limit_facets,
                              const Array<Element> & limit_segments,
                              const UInt & cent_node);

  /// use area as a weight
  Array<Element>
  findFacetsLoopByGraphByArea(const Array<Element> & limit_facets,
                              const Array<Element> & limit_segments,
                              const Array<Element> & preinserted_facets,
                              const UInt & cent_node);

  /// include only facets that have eff_stress >= 1
  Array<Element>
  findStressedFacetLoopAroundNode(const Array<Element> & limit_facets,
                                  const Array<Element> & limit_segments,
                                  const UInt & cent_node, Real min_dot);

  /// insert single facets before searching for the long loops
  Array<Element> findSingleFacetLoop(const Array<Element> & limit_facets,
                                     const Array<Element> & limit_segments,
                                     const UInt & cent_node);

  /// check if the cohesive can be inserted and nodes are not on the partition
  /// border and crack
  bool isFacetAndNodesGood(const Element & facet, UInt material_id = UInt(-1),
                           bool check_crack_facets = false);

  /// check if 2 facets are bounding a common solid element
  bool belong2SameElement(const Element & facet1, const Element & facet2);

  /// check if the node is surrounded by the material
  bool isNodeWithinMaterial(Element & node, UInt material_id);

  /// pick two neighbors of a central facet: returns true if success
  bool pickFacetNeighbors(Element & cent_facet);

  /// optimise doubled facets group, insert facets, and cohesive elements,
  /// update connectivities
  void insertOppositeFacetsAndCohesives();

  /// no cohesive elements in-between neighboring solid elements
  void preventCohesiveInsertionInNeighbors();

  /// assign crack tags to the clusters
  void assignCrackNumbers();

  /// change coordinates of central crack nodes to create an artificial gap
  void insertGap(const Real gap_ratio);

  /// same in 3D
  void insertGap3D(const Real gap_ratio);

  void onNodesAdded(const Array<UInt> & new_nodes,
                    const NewNodesEvent & nodes_event);

public:
  /// update the element group in case connectivity changed after cohesive
  /// elements insertion
  void updateElementGroup(const std::string group_name);

  /// insert cohesives on closed loops of stressed facets
  template <UInt dim> UInt insertLoopsOfStressedCohesives(Real min_dot);

  /// insert single cohesives based on effective stresses
  template <UInt dim> UInt insertStressedCohesives();

  /// find critical facets by looping on contour nodes
  std::map<UInt, std::map<UInt, UInt>> findCriticalFacetsOnContourByNodes(
      const std::map<Element, Element> & contour_subfacets_coh_el,
      const std::map<Element, UInt> & surface_subfacets_crack_nb,
      const std::map<UInt, std::set<Element>> & contour_nodes_subfacets,
      const std::set<UInt> & surface_nodes, Real av_stress_threshold);

  /// find critical facets by looping on contour nodes
  std::map<UInt, std::map<UInt, UInt>> findStressedFacetsLoops(
      const std::map<Element, Element> & contour_subfacets_coh_el,
      const std::map<Element, UInt> & surface_subfacets_crack_nb,
      const std::map<UInt, std::set<Element>> & contour_nodes_subfacets,
      const std::set<UInt> & surface_nodes, Real min_dot);

  /// average integrated eff stress over the area of the loop
  Real averageEffectiveStressInMultipleFacets(Array<Element> & facet_loop);

  void decomposeFacetLoopPerMaterial(
      const Array<Element> & facet_loop, UInt crack_nb,
      std::map<UInt, std::map<UInt, UInt>> & mat_index_facet_nbs_crack_nbs);

  /// apply eigen opening at all cohesives (including ghosts)
  void applyEigenOpening(Real eigen_strain);

  /// distribute total fluid volume along existing crack surface
  void applyEigenOpeningFluidVolumeBased(Real fluid_volume_ratio);

  /// apply pairs of equal opposite forces at the central nodes of pre-cracks
  void applyPointForceToCrackCentralNodes(Real force_norm);

  /// apply forces as before but with certain delay defined per crack
  void applyPointForceDelayed(Real loading_rate,
                              const Array<Real> loading_times, Real time,
                              Real multiplier = 1.);
  /// apply forces within an expanding sphere
  void applyPointForceDistributed(Real radius, Real F);

  /// apply force by equally splitting between first loop nodes
  void applyPointForcesFacetLoop(Real load);

  /// outputs crack area, volume into a file
  void outputCrackData(std::ofstream & file_output, Real time);

  /// output individual crack volumes
  void outputCrackVolumes(std::ofstream & file_output, Real time);

  /// computes crack area and volume per material
  std::tuple<Real, Real> computeCrackData(const ID & material_name);

  /// computes normal openings between crack central nodes and outputs them
  void outputCrackOpenings(std::ofstream & file_output, Real time);

  /// computes normal openings between crack central nodes and outputs them per
  /// processor
  void outputCrackOpeningsPerProc(std::ofstream & file_output, Real time);

  /// splits crack central nodes into pairs and computes averaged normals
  void identifyCrackCentralNodePairsAndNormals();

  /// compute crack contour segments, surface segments, contour nodes and
  /// surface nodes
  std::tuple<std::map<Element, Element>, std::map<Element, UInt>,
             std::set<UInt>, std::map<UInt, std::set<Element>>>
  determineCrackSurface();

  /// search for crack contour,etc. in a single facet
  void searchCrackSurfaceInSingleFacet(
      Element & facet, std::map<Element, Element> & contour_subfacets_coh_el,
      std::map<Element, UInt> & surface_subfacets_crack_nb,
      std::set<UInt> & surface_nodes,
      std::map<UInt, std::set<Element>> & contour_nodes_subfacets,
      std::set<Element> & visited_subfacets);

  /// update delta_max values in all sequential cohesive mats
  template <UInt dim> UInt updateDeltaMax();

  /// max ratio traction norm over max possible traction (SLA)
  template <UInt dim> std::tuple<Real, Element, UInt> getMaxDeltaMaxExcess();

  /// apply delta u on nodes
  void applyDeltaU(Real delta_u);

  /// apply eigenstrain on a specific material
  void applyEigenStrain(const Matrix<Real> & prestrain,
                        const ID & material_name = "gel");

  /// sets all facets within specified material to not be considered for
  /// insertion
  Real preventCohesiveInsertionInMaterial(std::string facet_mat_name);

  /// set update stiffness in all linear sequential cohesives
  template <UInt dim> void setUpdateStiffness(bool update_stiffness);
  /* ----------------------------------------------------------------- */

  /// RVE part

  /// apply boundary contions based on macroscopic deformation gradient
  virtual void
  applyBoundaryConditionsRve(const Matrix<Real> & displacement_gradient);

  /// apply homogeneous temperature field from the macroscale level to the
  /// RVEs
  virtual void applyHomogeneousTemperature(const Real & temperature);

  /// remove temperature from RVE
  virtual void removeTemperature();

  /// averages scalar field over the WHOLE!!! volume of a model
  Real averageScalarField(const ID & field_name);

  /// compute average stress or strain in the model
  Real averageTensorField(UInt row_index, UInt col_index,
                          const ID & field_type);

  /// compute effective stiffness of the RVE (in tension by default)
  void homogenizeStiffness(Matrix<Real> & C_macro, bool tensile_test = true);

  /// compute average eigenstrain
  void homogenizeEigenGradU(Matrix<Real> & eigen_gradu_macro);

  /// compute damage to the RVE area ratio
  void computeDamageRatio(Real & damage);

  /// compute damage to material area ratio
  void computeDamageRatioPerMaterial(Real & damage_ratio,
                                     const ID & material_name);

  /// compute ratio between total crack and RVE volumes
  void computeCrackVolume(Real & crack_volume_ratio);

  /// compute crack volume to material volume ratio
  void computeCrackVolumePerMaterial(Real & crack_volume,
                                     const ID & material_name);

  /// dump the RVE
  void dumpRve();

  /// dump the RVE with the dump number
  void dumpRve(UInt dump_nb);

  /// compute average stress in the RVE
  void homogenizeStressField(Matrix<Real> & stress);

  /// compute hydrostatic part of the stress and assign homogen direction
  void setStiffHomogenDir(Matrix<Real> & stress);

  /// storing nodal fields before tests
  void storeNodalFields();

  /// restoring nodal fields before tests
  void restoreNodalFields();

  /// resetting nodal fields
  void resetNodalFields();

  /// restoring internal fields after tests
  void restoreInternalFields();

  /// resetting internal fields with history
  void resetInternalFields();

  /// store damage in materials who have damage
  void storeDamageField();

  /// assign stored values to the current ones
  void restoreDamageField();

  /// find the corner nodes
  void findCornerNodes();

  /// perform virtual testing
  void performVirtualTesting(const Matrix<Real> & H,
                             Matrix<Real> & eff_stresses,
                             Matrix<Real> & eff_strains, const UInt test_no);

  /// clear the eigenstrain (to exclude stresses due to internal pressure)
  void clearEigenStrain(const ID & material_name = "gel");

  /// store elemental eigenstrain in an array
  void storeEigenStrain(Array<Real> & stored_eig,
                        const ID & material_name = "gel");

  /// restore eigenstrain  from previously stored values
  void restoreEigenStrain(Array<Real> & stored_eig,
                          const ID & material_name = "gel");
  /* ----------------------------------------------------------------  */
public:
  // Accessors
  bool isTensileHomogen() { return this->tensile_homogenization; };

  /// phase volumes
  Real getPhaseVolume(const std::string & material_name) {
    if (not this->phase_volumes.size()) {
      computePhaseVolumes();
    }
    return this->phase_volumes.find(material_name)->second;
  };

  /// set the value of the insertion flag
  AKANTU_SET_MACRO(CohesiveInsertion, cohesive_insertion, bool);

  AKANTU_GET_MACRO(NodesEffStress, nodes_eff_stress, const Array<Real> &);
  /// get the corner nodes
  AKANTU_GET_MACRO(CornerNodes, corner_nodes, const Array<UInt> &);

  AKANTU_GET_MACRO(CrackFacets, crack_facets, const Array<Array<Element>> &);
  AKANTU_GET_MACRO(CrackFacetsFromMesh, crack_facets_from_mesh,
                   const Array<Array<Element>> &);
  AKANTU_GET_MACRO(CrackCentralNodes, crack_central_nodes, const Array<UInt> &);
  AKANTU_GET_MACRO(CrackCentralNodePairs, crack_central_node_pairs,
                   const auto &);
  AKANTU_GET_MACRO(CrackNormalsPairs, crack_normals_pairs, const auto &);

  /* --------------------------------------------------------------------- */
  /* Members */
  /* --------------------------------------------------------------------- */

private:
  SolidMechanicsModel & model;

protected:
  /// volume of the RVE
  Real volume;

  /// corner nodes 1, 2, 3, 4 (see Leonardo's thesis, page 98)
  Array<UInt> corner_nodes;

  /// bottom nodes
  std::unordered_set<UInt> bottom_nodes;

  /// left nodes
  std::unordered_set<UInt> left_nodes;

  /// dump counter
  UInt nb_dumps;

  /// flag to activate expansion through cohesive elements
  bool cohesive_insertion;

  /// booleans for applying delta u
  bool doubled_facets_ready;
  bool crack_central_nodes_ready;

  // arrays to store nodal values during virtual tests
  Array<Real> disp_stored;
  Array<Real> ext_force_stored;
  Array<bool> boun_stored;

  /// if stiffness homogenization will be done in tension
  bool tensile_homogenization{false};

  /// phase volumes
  std::map<std::string, Real> phase_volumes;

  /// array to store flags on nodes position modification
  Array<bool> modified_pos;

  /// array to store flags on nodes that are synchronized between processors
  Array<bool> partition_border_nodes;

  /// average eff stress on facet loops around a node
  Array<Real> nodes_eff_stress;

  /// array to store flags on nodes where cracks are inserted
  Array<bool> crack_nodes;

  /// Vector to store the initially inserted facets per single crack
  Array<Array<Element>> crack_facets;
  Array<Array<Element>> crack_facets_from_mesh;
  /// duplicated nodes. first 1/2 - initialnodes, second - duplicated
  Array<UInt> crack_central_nodes;
  /// crack central nodes devided into pairs
  Array<std::pair<UInt, UInt>> crack_central_node_pairs;
  /// normals per node averaged per crack of each site
  Array<std::pair<Vector<Real>, Vector<Real>>> crack_normals_pairs;
};

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
/* ------------------------------------------------------------------- */
/* Boundary conditions functors */
/* ------------------------------------------------------------------- */

class Pressure : public BC::Neumann::NeumannFunctor {
public:
  Pressure(SolidMechanicsModel & model, const Array<Real> & pressure_on_qpoint,
           const Array<Real> & quad_coords)
      : model(model), pressure_on_qpoint(pressure_on_qpoint),
        quad_coords(quad_coords) {}

  inline void operator()(const IntegrationPoint & quad_point,
                         Vector<Real> & dual, const Vector<Real> & coord,
                         const Vector<Real> & normal) const {

    // get element types
    auto && mesh = model.getMesh();
    const UInt dim = mesh.getSpatialDimension();
    const GhostType gt = akantu::_not_ghost;
    const UInt facet_nb = quad_point.element;
    const ElementKind cohesive_kind = akantu::_ek_cohesive;
    const ElementType type_facet = quad_point.type;
    const ElementType type_coh = FEEngine::getCohesiveElementType(type_facet);
    auto && cohesive_conn = mesh.getConnectivity(type_coh, gt);
    const UInt nb_nodes_coh_elem = cohesive_conn.getNbComponent();
    auto && facet_conn = mesh.getConnectivity(type_facet, gt);
    const UInt nb_nodes_facet = facet_conn.getNbComponent();
    auto && fem_boundary = model.getFEEngineBoundary();
    UInt nb_quad_points = fem_boundary.getNbIntegrationPoints(type_facet, gt);
    auto facet_nodes_it = make_view(facet_conn, nb_nodes_facet).begin();

    AKANTU_DEBUG_ASSERT(nb_nodes_coh_elem == 2 * nb_nodes_facet,
                        "Different number of nodes belonging to one cohesive "
                        "element face and facet");
    // const akantu::Mesh &mesh_facets = cohesive_mesh.getMeshFacets();
    const Array<std::vector<Element>> & elem_to_subelem =
        mesh.getElementToSubelement(type_facet, gt);
    const auto & adjacent_elems = elem_to_subelem(facet_nb);
    auto normal_corrected = normal;

    // loop over all adjacent elements
    UInt coh_elem_nb;
    if (not adjacent_elems.empty()) {
      for (UInt f : arange(adjacent_elems.size())) {
        if (adjacent_elems[f].kind() != cohesive_kind)
          continue;
        coh_elem_nb = adjacent_elems[f].element;
        Array<UInt> upper_nodes(nb_nodes_coh_elem / 2);
        Array<UInt> lower_nodes(nb_nodes_coh_elem / 2);
        for (UInt node : arange(nb_nodes_coh_elem / 2)) {
          upper_nodes(node) = cohesive_conn(coh_elem_nb, node);
          lower_nodes(node) =
              cohesive_conn(coh_elem_nb, node + nb_nodes_coh_elem / 2);
        }
        bool upper_face = true;
        bool lower_face = true;
        Vector<UInt> facet_nodes = facet_nodes_it[facet_nb];
        for (UInt s : arange(nb_nodes_facet)) {
          auto idu = upper_nodes.find(facet_nodes(s));
          auto idl = lower_nodes.find(facet_nodes(s));
          if (idu == UInt(-1))
            upper_face = false;
          else if (idl == UInt(-1))
            lower_face = false;
        }
        if (upper_face && not lower_face)
          normal_corrected *= -1;
        else if (not upper_face && lower_face)
          normal_corrected *= 1;
        else
          AKANTU_EXCEPTION("Error in defining side of the cohesive element");
        break;
      }
    }

    auto flow_qpoint_it = make_view(quad_coords, dim).begin();
    bool node_found;
    for (auto qp : arange(nb_quad_points)) {
      const Vector<Real> flow_quad_coord =
          flow_qpoint_it[coh_elem_nb * nb_quad_points + qp];
      if (flow_quad_coord != coord)
        continue;
      Real P = pressure_on_qpoint(coh_elem_nb * nb_quad_points + qp);
      dual = P * normal_corrected;
      node_found = true;
    }
    if (not node_found)
      AKANTU_EXCEPTION("Quad point is not found in the flow mesh");
  }

protected:
  SolidMechanicsModel & model;
  const Array<Real> & pressure_on_qpoint;
  const Array<Real> & quad_coords;
};

/* ------------------------------------------------------------------ */
class PressureSimple : public BC::Neumann::NeumannFunctor {
public:
  PressureSimple(SolidMechanicsModel & model, const Real pressure,
                 const std::string group_name)
      : model(model), pressure(pressure), group_name(group_name) {}

  inline void operator()(const IntegrationPoint & quad_point,
                         Vector<Real> & dual, const Vector<Real> & /*coord*/,
                         const Vector<Real> & normal) const {

    // get element types
    auto && mesh = model.getMesh();
    AKANTU_DEBUG_ASSERT(mesh.elementGroupExists(group_name),
                        "Element group is not registered in the mesh");
    const GhostType gt = akantu::_not_ghost;
    const UInt facet_nb = quad_point.element;
    const ElementType type_facet = quad_point.type;
    auto && facet_conn = mesh.getConnectivity(type_facet, gt);
    const UInt nb_nodes_facet = facet_conn.getNbComponent();
    auto facet_nodes_it = make_view(facet_conn, nb_nodes_facet).begin();
    auto & group = mesh.getElementGroup(group_name);
    Array<UInt> element_ids = group.getElements(type_facet);
    AKANTU_DEBUG_ASSERT(element_ids.size(),
                        "Provided group doesn't contain this element type");
    auto id = element_ids.find(facet_nb);
    AKANTU_DEBUG_ASSERT(id != UInt(-1),
                        "Quad point doesn't belong to this element group");

    auto normal_corrected = normal;

    if (id < element_ids.size() / 2)
      normal_corrected *= -1;
    else if (id >= element_ids.size())
      AKANTU_EXCEPTION("Error in defining side of the cohesive element");
    else
      normal_corrected *= 1;

    dual = pressure * normal_corrected;
  }

protected:
  SolidMechanicsModel & model;
  const Real pressure;
  const std::string group_name;
};

/* ------------------------------------------------------------------- */
class PressureVolumeDependent : public BC::Neumann::NeumannFunctor {
public:
  PressureVolumeDependent(SolidMechanicsModel & model,
                          const Real fluid_volume_ratio,
                          const std::string group_name,
                          const Real compressibility)
      : model(model), fluid_volume_ratio(fluid_volume_ratio),
        group_name(group_name), compressibility(compressibility) {}

  inline void operator()(const IntegrationPoint & quad_point,
                         Vector<Real> & dual, const Vector<Real> & /*coord*/,
                         const Vector<Real> & /*normal*/) const {

    // get element types
    auto && mesh = model.getMesh();
    AKANTU_DEBUG_ASSERT(mesh.elementGroupExists(group_name),
                        "Element group is not registered in the mesh");
    auto dim = mesh.getSpatialDimension();
    const GhostType gt = akantu::_not_ghost;
    const UInt facet_nb = quad_point.element;
    const ElementType type_facet = quad_point.type;
    auto && facet_conn = mesh.getConnectivity(type_facet, gt);
    const UInt nb_nodes_facet = facet_conn.getNbComponent();
    auto facet_nodes_it = make_view(facet_conn, nb_nodes_facet).begin();
    auto & group = mesh.getElementGroup(group_name);
    Array<UInt> element_ids = group.getElements(type_facet);
    auto && pos = mesh.getNodes();
    const auto pos_it = make_view(pos, dim).begin();
    auto && disp = model.getDisplacement();
    const auto disp_it = make_view(disp, dim).begin();
    auto && fem_boundary = model.getFEEngineBoundary();
    UInt nb_quad_points = fem_boundary.getNbIntegrationPoints(type_facet, gt);

    AKANTU_DEBUG_ASSERT(element_ids.size(),
                        "Provided group doesn't contain this element type");
    auto id = element_ids.find(facet_nb);
    AKANTU_DEBUG_ASSERT(id != UInt(-1),
                        "Quad point doesn't belong to this element group");

    // get normal to the current positions
    const auto & current_pos = model.getCurrentPosition();
    Array<Real> quad_normals(0, dim);
    fem_boundary.computeNormalsOnIntegrationPoints(current_pos, quad_normals,
                                                   type_facet, gt);
    auto normals_it = quad_normals.begin(dim);
    Vector<Real> normal_corrected(
        normals_it[quad_point.element * nb_quad_points + quad_point.num_point]);

    // auto normal_corrected = normal;
    UInt opposite_facet_nb(-1);
    if (id < element_ids.size() / 2) {
      normal_corrected *= -1;
      opposite_facet_nb = element_ids(id + element_ids.size() / 2);
    } else if (id >= element_ids.size())
      AKANTU_EXCEPTION("Error in defining side of the cohesive element");
    else {
      normal_corrected *= 1;
      opposite_facet_nb = element_ids(id - element_ids.size() / 2);
    }

    /// compute current area of a gap
    Vector<UInt> first_facet_nodes = facet_nodes_it[facet_nb];
    Vector<UInt> second_facet_nodes = facet_nodes_it[opposite_facet_nb];
    /* corners of a quadrangle consequently
            A---M---B
            \      |
            D--N--C */
    UInt A, B, C, D;
    A = first_facet_nodes(0);
    B = first_facet_nodes(1);
    C = second_facet_nodes(0);
    D = second_facet_nodes(1);

    /// quadrangle's area through diagonals
    Vector<Real> AC, BD;
    Vector<Real> A_pos = Vector<Real>(pos_it[A]) + Vector<Real>(disp_it[A]);
    Vector<Real> B_pos = Vector<Real>(pos_it[B]) + Vector<Real>(disp_it[B]);
    Vector<Real> C_pos = Vector<Real>(pos_it[C]) + Vector<Real>(disp_it[C]);
    Vector<Real> D_pos = Vector<Real>(pos_it[D]) + Vector<Real>(disp_it[D]);
    Vector<Real> M_pos = (A_pos + B_pos) * 0.5;
    Vector<Real> N_pos = (C_pos + D_pos) * 0.5;
    Vector<Real> MN = M_pos - N_pos;
    Vector<Real> AB = A_pos - B_pos;
    Vector<Real> AB_0 = Vector<Real>(pos_it[A]) - Vector<Real>(pos_it[B]);

    // fluid volume computed as AB * thickness (AB * ratio)
    Real fluid_volume = AB_0.norm() * AB_0.norm() * fluid_volume_ratio;
    Real current_volume = AB.norm() * MN.norm();
    current_volume =
        Math::are_float_equal(current_volume, 0.) ? 0 : current_volume;
    Real volume_change = current_volume - fluid_volume;
    Real pressure_change{0};
    if (volume_change < 0) {
      pressure_change = -volume_change / fluid_volume / this->compressibility;
    }
    dual = pressure_change * normal_corrected;
  }

protected:
  SolidMechanicsModel & model;
  const Real fluid_volume_ratio;
  const std::string group_name;
  const Real compressibility;
};

/* ------------------------------------------------------------------- */
class PressureVolumeDependent3D : public BC::Neumann::NeumannFunctor {
public:
  PressureVolumeDependent3D(SolidMechanicsModelCohesive & model,
                            const Real fluid_volume_ratio,
                            const Real compressibility,
                            const Array<Array<Element>> crack_facets_from_mesh)
      : model(model), fluid_volume_ratio(fluid_volume_ratio),
        compressibility(compressibility),
        crack_facets_from_mesh(crack_facets_from_mesh) {}

  inline void operator()(const IntegrationPoint & quad_point,
                         Vector<Real> & dual, const Vector<Real> & /*coord*/,
                         const Vector<Real> & /*normal*/) const {

    // get element types
    auto && mesh = model.getMesh();
    const FEEngine & fe_engine = model.getFEEngine("CohesiveFEEngine");
    auto dim = mesh.getSpatialDimension();
    Element facet{quad_point.type, quad_point.element, quad_point.ghost_type};
    ElementType type_cohesive = FEEngine::getCohesiveElementType(facet.type);
    auto nb_quad_coh = fe_engine.getNbIntegrationPoints(type_cohesive);
    auto && facet_nodes = mesh.getConnectivity(facet);
    auto && fem_boundary = model.getFEEngineBoundary();
    UInt nb_quad_points =
        fem_boundary.getNbIntegrationPoints(facet.type, facet.ghost_type);

    // get normal to the current positions
    const auto & current_pos = model.getCurrentPosition();
    Array<Real> quad_normals(0, dim);
    fem_boundary.computeNormalsOnIntegrationPoints(
        current_pos, quad_normals, facet.type, facet.ghost_type);
    auto normals_it = quad_normals.begin(dim);
    Vector<Real> normal(
        normals_it[quad_point.element * nb_quad_points + quad_point.num_point]);
    // search for this facet in the crack facets array
    UInt loading_site_nb(-1);
    UInt id_in_array(-1);
    for (auto && data : enumerate(this->crack_facets_from_mesh)) {
      auto & one_site_facets = std::get<1>(data);
      id_in_array = one_site_facets.find(facet);
      if (id_in_array != UInt(-1)) {
        loading_site_nb = std::get<0>(data);
        break;
      }
    }
    AKANTU_DEBUG_ASSERT(loading_site_nb != UInt(-1),
                        "Quad point doesn't belong to the loading facets");

    // find connected solid element
    auto && facet_neighbors = mesh.getElementToSubelement(facet);
    Element solid_el{ElementNull};
    for (auto && neighbor : facet_neighbors) {
      if (mesh.getKind(neighbor.type) == _ek_regular) {
        solid_el = neighbor;
        break;
      }
    }
    AKANTU_DEBUG_ASSERT(solid_el != ElementNull,
                        "Couldn't identify neighboring solid el");

    // find the out-of-plane solid node
    auto && solid_nodes = mesh.getConnectivity(solid_el);

    UInt out_node(-1);
    for (auto && solid_node : solid_nodes) {
      auto ret = std::find(facet_nodes.begin(), facet_nodes.end(), solid_node);
      if (ret == facet_nodes.end()) {
        out_node = solid_node;
        break;
      }
    }
    AKANTU_DEBUG_ASSERT(out_node != UInt(-1),
                        "Couldn't identify out-of-plane node");

    // build the reference vector
    Vector<Real> ref_vector(dim);
    mesh.getBarycenter(facet, ref_vector);
    Vector<Real> pos(mesh.getNodes().begin(dim)[out_node]);
    ref_vector = pos - ref_vector;

    // check if ref vector and the normal are in the same half space
    if (ref_vector.dot(normal) < 0) {
      normal *= -1;
    }

    // compute volume (area * normal_opening) of current fluid body
    // form cohesive element filter from the 1st half of facet filter
    std::set<Element> site_cohesives;
    for (auto & crack_facet : this->crack_facets_from_mesh(loading_site_nb)) {
      if (crack_facet.type == facet.type and
          crack_facet.ghost_type == facet.ghost_type) {
        // find connected cohesive
        auto & connected_els = mesh.getElementToSubelement(crack_facet);
        for (auto & connected_el : connected_els) {
          if (connected_el.type == type_cohesive) {
            site_cohesives.emplace(connected_el);
            break;
          }
        }
      }
    }
    // integrate normal opening over identified element filter
    Real site_volume{0};
    Real site_area{0};
    const Array<UInt> & material_index_vec =
        model.getMaterialByElement(type_cohesive, facet.ghost_type);
    const Array<UInt> & material_local_numbering_vec =
        model.getMaterialLocalNumbering(type_cohesive, facet.ghost_type);
    for (auto & coh_el : site_cohesives) {
      Material & material =
          model.getMaterial(material_index_vec(coh_el.element));
      UInt material_local_num = material_local_numbering_vec(coh_el.element);
      Array<UInt> single_el_array;
      single_el_array.push_back(coh_el.element);
      auto & opening_norm_array = material.getInternal<Real>(
          "normal_opening_norm")(coh_el.type, coh_el.ghost_type);
      Array<Real> opening_norm_el;
      for (UInt i = 0; i != nb_quad_coh; i++) {
        auto & opening_per_quad =
            opening_norm_array(material_local_num * nb_quad_coh + i);
        opening_norm_el.push_back(opening_per_quad);
      }

      site_volume += fe_engine.integrate(opening_norm_el, coh_el.type,
                                         coh_el.ghost_type, single_el_array);
      Array<Real> area(fe_engine.getNbIntegrationPoints(type_cohesive), 1, 1.);
      site_area += fe_engine.integrate(area, coh_el.type, coh_el.ghost_type,
                                       single_el_array);
    }

    Real fluid_volume = site_area * fluid_volume_ratio;
    Real volume_change = site_volume - fluid_volume;

    Real pressure_change{0};
    if (volume_change < 0) {
      pressure_change = -volume_change / fluid_volume / this->compressibility;
    }
    std::cout << " volume = " << site_area << " x " << fluid_volume_ratio
              << " = " << fluid_volume << " site volume " << site_volume
              << " volume change " << volume_change << " pressure delta "
              << pressure_change << std::endl;
    dual = pressure_change * normal;
  }

protected:
  SolidMechanicsModelCohesive & model;
  const Real fluid_volume_ratio;
  const Real compressibility;
  const Array<Array<Element>> crack_facets_from_mesh;
};

/* ------------------------------------------------------------------ */
class DeltaU : public BC::Dirichlet::DirichletFunctor {
public:
  DeltaU(const SolidMechanicsModel & model, const Real delta_u,
         const Array<std::tuple<UInt, UInt>> & node_pairs)
      : model(model), delta_u(delta_u), node_pairs(node_pairs),
        displacement(model.getDisplacement()) {}

  inline void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal,
                         const Vector<Real> &) const {

    // get element types
    auto && mesh = model.getMesh();
    const UInt dim = mesh.getSpatialDimension();
    auto && mesh_facets = mesh.getMeshFacets();
    auto disp_it = make_view(displacement, dim).begin();
    CSR<Element> nodes_to_elements;
    MeshUtils::buildNode2Elements(mesh_facets, nodes_to_elements, dim - 1);

    // get actual distance between two nodes
    Vector<Real> node_disp(disp_it[node]);
    Vector<Real> other_node_disp(dim);
    bool upper_face = false;
    bool lower_face = false;
    for (auto && pair : this->node_pairs) {
      if (node == std::get<0>(pair)) {
        other_node_disp = disp_it[std::get<1>(pair)];
        upper_face = true;
        break;
      } else if (node == std::get<1>(pair)) {
        other_node_disp = disp_it[std::get<0>(pair)];
        lower_face = true;
        break;
      }
    }
    AKANTU_DEBUG_ASSERT(upper_face == true or lower_face == true,
                        "Error in identifying the node in tuple");
    Real sign = -upper_face + lower_face;

    // compute normal at node (averaged between two surfaces)
    Vector<Real> normal(dim);
    for (auto & elem : nodes_to_elements.getRow(node)) {
      if (mesh.getKind(elem.type) != _ek_regular)
        continue;
      if (elem.ghost_type != _not_ghost)
        continue;
      auto & doubled_facets_array = mesh_facets.getData<bool>(
          "doubled_facets", elem.type, elem.ghost_type);
      if (doubled_facets_array(elem.element) != true)
        continue;

      auto && fe_engine_facet = model.getFEEngine("FacetsFEEngine");
      auto nb_qpoints_per_facet =
          fe_engine_facet.getNbIntegrationPoints(elem.type, elem.ghost_type);
      const auto & normals_on_quad =
          fe_engine_facet.getNormalsOnIntegrationPoints(elem.type,
                                                        elem.ghost_type);
      auto normals_it = make_view(normals_on_quad, dim).begin();
      normal +=
          sign * Vector<Real>(normals_it[elem.element * nb_qpoints_per_facet]);
    }
    normal /= normal.norm();

    // get distance between two nodes in normal direction
    Real node_disp_norm = node_disp.dot(normal);
    Real other_node_disp_norm = other_node_disp.dot(-1. * normal);
    Real dist = node_disp_norm + other_node_disp_norm;
    Real prop_factor = dist == 0. ? 0.5 : node_disp_norm / dist;

    // get correction displacement
    Real correction = delta_u - dist;

    // apply absolute value of displacement
    primal += normal * correction * prop_factor;
    flags.set(false);
  }

protected:
  const SolidMechanicsModel & model;
  const Real delta_u;
  const Array<std::tuple<UInt, UInt>> node_pairs;
  Array<Real> displacement;
};

/* ------------------------------------------------------------------ */
/* solid mechanics model cohesive + delta d at the crack nodes */
/* ------------------------------------------------------------------ */
class SolidMechanicsModelCohesiveDelta : public SolidMechanicsModelCohesive {
public:
  SolidMechanicsModelCohesiveDelta(Mesh & mesh)
      : SolidMechanicsModelCohesive(mesh), mesh(mesh), delta_u(0.) {}

  /* ------------------------------------------------------------------*/
  void assembleInternalForces() {
    // displacement correction
    applyDisplacementDifference();

    SolidMechanicsModelCohesive::assembleInternalForces();
  }
  /* ----------------------------------------------------------------- */
  void applyDisplacementDifference() {
    auto dim = this->mesh.getSpatialDimension();
    auto & disp = this->getDisplacement();
    auto & boun = this->getBlockedDOFs();

    // get normal to the initial positions
    auto it_disp = make_view(disp, dim).begin();
    auto it_boun = make_view(boun, dim).begin();

    for (auto && data : zip(crack_central_node_pairs, crack_normals_pairs)) {
      auto && node_pair = std::get<0>(data);
      auto && normals_pair = std::get<1>(data);

      auto node1 = node_pair.first;
      auto node2 = node_pair.second;
      auto normal1 = normals_pair.first;
      auto normal2 = normals_pair.second;
      Vector<Real> displ1 = it_disp[node1];
      Vector<Real> displ2 = it_disp[node2];
      if (mesh.isPeriodicSlave(node1)) {
        displ1.copy(displ2);
      } else {
        displ2.copy(displ1);
      }
      displ1 += normal1 * this->delta_u / 2.;
      displ2 += normal2 * this->delta_u / 2.;
      // Vector<bool> dof = it_boun[node];
      // dof.set(true);
    }
  }
  /* ------------------------------------------------------------------*/
public:
  // Acessors
  AKANTU_GET_MACRO_NOT_CONST(DeltaU, delta_u, Real &);
  using NodePairsArray = Array<std::pair<UInt, UInt>>;
  AKANTU_SET_MACRO(CrackNodePairs, crack_central_node_pairs, NodePairsArray);
  using NormalPairsArray = Array<std::pair<Vector<Real>, Vector<Real>>>;
  AKANTU_SET_MACRO(CrackNormalsPairs, crack_normals_pairs, NormalPairsArray);

protected:
  Mesh & mesh;
  Real delta_u;
  Array<std::pair<UInt, UInt>> crack_central_node_pairs;
  Array<std::pair<Vector<Real>, Vector<Real>>> crack_normals_pairs;
};
/* ------------------------------------------------------------------ */

} // namespace akantu

#endif /* __AKANTU_RVE_TOOLS_HH__ */
