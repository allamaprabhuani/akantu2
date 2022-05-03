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

} // namespace akantu

#endif /* __AKANTU_RVE_TOOLS_HH__ */
