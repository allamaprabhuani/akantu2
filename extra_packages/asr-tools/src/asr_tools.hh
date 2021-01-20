/**
 * @file   asr_tools.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue Jan 16 10:26:53 2014
 *
 * @brief  tools for the analysis of ASR samples
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
#ifdef AKANTU_COHESIVE_ELEMENT
#include "solid_mechanics_model_cohesive.hh"
#endif
#ifdef AKANTU_FLUID_DIFFUSION
#include "fluid_diffusion_model.hh"
#endif
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ASR_TOOLS_HH__
#define __AKANTU_ASR_TOOLS_HH__

namespace akantu {
class NodeGroup;
class SolidMechanicsModel;
} // namespace akantu

namespace akantu {

class ASRTools : public MeshEventHandler {
public:
  ASRTools(SolidMechanicsModel & model);

  virtual ~ASRTools() = default;

  /// This function is used in case of uni-axial boundary conditions:
  /// It detects the nodes, on which a traction needs to be applied and stores
  /// them in a node group
  void fillNodeGroup(NodeGroup & node_group, bool multi_axial = false);

  /// compute volumes of each phase
  void computePhaseVolumes();

  /// compute volume of the model passed to ASRTools (RVE if FE2_mat)
  void computeModelVolume();

  /// Apply free expansion boundary conditions
  void applyFreeExpansionBC();

  /// Apply loaded boundary conditions
  void applyLoadedBC(const Vector<Real> & traction, const ID & element_group,
                     bool multi_axial);

  /// This function computes the average displacement along one side of the
  /// sample,
  /// caused by the expansion of gel pockets
  Real computeAverageDisplacement(SpatialDirection direction);

  /// This function computes a compoment of average volumetric strain tensor,
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
                                             Real time);

  /// Averages strains & displacements + damage ratios in FE2 material
  void computeAveragePropertiesFe2Material(std::ofstream & file_output,
                                           Real time);

  /// Save the state of the simulation for subsequent restart
  void saveState(UInt current_step);

  /// Load to previous state of the simulation for restart
  bool loadState(UInt & current_step);

  /// compute the volume occupied by a given material
  Real computePhaseVolume(const ID & mat_name);

  /// apply eigenstrain in cracks, that are filled with gel
  void applyEigenGradUinCracks(const Matrix<Real> prescribed_eigen_grad_u,
                               const ID & material_name);

  /// apply eigenstrain in the cracks, that advanced in the last step
  void applyEigenGradUinCracks(const Matrix<Real> prescribed_eigen_grad_u,
                               const ElementTypeMapUInt & critical_elements,
                               bool to_add = false);

  /// fill fully cracked elements with gel
  void fillCracks(ElementTypeMapReal & saved_damage);

  /// drain the gel in fully cracked elements
  void drainCracks(const ElementTypeMapReal & saved_damage);

  Real computeSmallestElementSize();

  // /// apply homogeneous temperature on Solid Mechanics Model
  // void applyTemperatureFieldToSolidmechanicsModel(const Real & temperature);

  /// compute ASR strain by a sigmoidal rule (Larive, 1998)
  void computeASRStrainLarive(const Real & delta_time_day, const Real & T,
                              Real & ASRStrain, const Real & eps_inf,
                              const Real & time_ch_ref,
                              const Real & time_lat_ref, const Real & U_C,
                              const Real & U_L, const Real & T_ref);

  /// compute increase in gel strain within 1 timestep
  Real computeDeltaGelStrainThermal(const Real delta_time_day, const Real k,
                                    const Real activ_energy, const Real R,
                                    const Real T,
                                    Real & amount_reactive_particles,
                                    const Real saturation_const);

  /// compute linear increase in gel strain
  Real computeDeltaGelStrainLinear(const Real delta_time, const Real k);

  /// insert multiple blocks of cohesive elements
  void insertASRCohesivesRandomly(const UInt & nb_coh_elem,
                                  std::string matrix_mat_name, Real gap_ratio);
  /// insert multiple blocks of cohesive elements
  void insertASRCohesivesRandomly3D(const UInt & nb_coh_elem,
                                    std::string matrix_mat_name,
                                    Real gap_ratio);

  /// ASR insertion for 1st order 3D elements
  void insertASRCohesiveLoops3D(const UInt & nb_insertions,
                                std::string facet_mat_name, Real gap_ratio);

  /// insert block of cohesive elements based on the coord of the central
  void insertASRCohesivesByCoords(const Matrix<Real> & positions,
                                  Real gap_ratio = 0);

  /// communicates crack numbers from not ghosts to ghosts cohesive elements
  void communicateCrackNumbers();

protected:
  /// put flags on all nodes who have a ghost counterpart
  void communicateFlagsOnNodes();

  /// pick facets by passed coordinates
  void pickFacetsByCoord(const Matrix<Real> & positions);

  /// pick facets randomly within specified material
  void pickFacetsRandomly(UInt nb_insertions, std::string facet_mat_name);

  /// build closed facets loop around random point of specified material and
  /// returns number of successfully inserted ASR sites
  UInt closedFacetsLoopAroundPoint(UInt nb_insertions, std::string mat_name);

  /// check if the cohesive can be inserted and nodes are not on the partition
  /// border and ASR nodes
  bool isFacetAndNodesGood(const Element & facet, UInt material_id);

  /// check if the node is surrounded by the material
  bool isNodeWithinMaterial(Element & node, UInt material_id);

  /// pick two neighbors of a central facet in 2D: returns true if success
  bool pickFacetNeighborsOld(Element & cent_facet);

  /// version working for both 2d and 3d
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

  /// on elements added for asr-tools
  void onElementsAdded(const Array<Element> & elements,
                       const NewElementsEvent & element_event);

  void onNodesAdded(const Array<UInt> & new_nodes,
                    const NewNodesEvent & nodes_event);

public:
  /// update the element group in case connectivity changed after cohesive
  /// elements insertion
  void updateElementGroup(const std::string group_name);

  /// works only for the MaterialCohesiveLinearSequential
  template <UInt dim> UInt insertCohesiveElementsSelectively();

  /// insert multiple cohesives on contour
  template <UInt dim> UInt insertCohesiveElementsOnContour();

  /// apply eigen opening at all cohesives (including ghosts)
  template <UInt dim> void applyEigenOpening(Real eigen_strain);

  /// outputs crack area, volume into a file
  void outputCrackData(std::ofstream & file_output, Real time);

  /// computes crack area and volume per material
  std::tuple<Real, Real> computeCrackData(const ID & material_name);

  // /// apply self-weight force
  // void applyBodyForce();

  /// apply delta u on nodes
  void applyDeltaU(Real delta_u);

  /// apply eigenstrain on gel material
  void applyGelStrain(const Matrix<Real> & prestrain);
  /* ------------------------------------------------------------------------
   */
  /// RVE part

  /// apply boundary contions based on macroscopic deformation gradient
  virtual void
  applyBoundaryConditionsRve(const Matrix<Real> & displacement_gradient);

  /// apply homogeneous temperature field from the macroscale level to the
  /// RVEs
  virtual void applyHomogeneousTemperature(const Real & temperature);

  /// remove temperature from RVE on the end of ASR advancement
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
  void clearASREigenStrain();

  /// store elemental eigenstrain in an array
  void storeASREigenStrain(Array<Real> & stored_eig);

  /// restore eigenstrain in ASR sites from previously stored values
  void restoreASREigenStrain(Array<Real> & stored_eig);
  /* ------------------------------------------------------------------------
   */
public:
  // Accessors
  bool isTensileHomogen() { return this->tensile_homogenization; };

  /// phase volumes
  Real getPhaseVolume(const std::string & material_name) {
    if (not this->phase_volumes.size())
      computePhaseVolumes();

    return this->phase_volumes.find(material_name)->second;
  };

  /// set the value of the insertion flag
  AKANTU_SET_MACRO(CohesiveInsertion, cohesive_insertion, bool);

  /// get the corner nodes
  AKANTU_GET_MACRO(CornerNodes, corner_nodes, const Array<UInt> &);

  /* --------------------------------------------------------------------- */
  /* Members */
  /* --------------------------------------------------------------------- */

private:
  SolidMechanicsModel & model;

  // /// 2D hardcoded - no 3D support currently
  // using voigt_h = VoigtHelper<2>;

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

  /// flag to activate ASR expansion through cohesive elements
  bool cohesive_insertion;

  /// booleans for applying delta u
  bool doubled_facets_ready;
  bool doubled_nodes_ready;

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

  /// array to store flags on nodes where ASR elements are inserted
  Array<bool> ASR_nodes;
};

/* --------------------------------------------------------------------------
 */
/* ASR material selector */
/* --------------------------------------------------------------------------
 */
class GelMaterialSelector : public MeshDataMaterialSelector<std::string> {
public:
  GelMaterialSelector(SolidMechanicsModel & model, std::string gel_material,
                      const UInt nb_gel_pockets,
                      std::string aggregate_material = "aggregate",
                      bool gel_pairs = false)
      : MeshDataMaterialSelector<std::string>("physical_names", model),
        model(model), gel_material(gel_material),
        nb_gel_pockets(nb_gel_pockets), nb_placed_gel_pockets(0),
        aggregate_material(aggregate_material), gel_pairs(gel_pairs) {}

  void initGelPocket() {
    aggregate_material_id = model.getMaterialIndex(aggregate_material);

    Mesh & mesh = this->model.getMesh();
    UInt dim = model.getSpatialDimension();
    //    Element el{_triangle_3, 0, _not_ghost};
    for (auto el_type : model.getMaterial(aggregate_material)
                            .getElementFilter()
                            .elementTypes(dim)) {

      const auto & filter =
          model.getMaterial(aggregate_material).getElementFilter()(el_type);
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

      /// generate the gel pockets
      srand(0.);
      Vector<Real> center(dim);
      std::set<int> checked_baries;
      while (nb_placed_gel_pockets != nb_gel_pockets) {
        /// get a random bary center
        UInt bary_id = rand() % nb_element;
        if (checked_baries.find(bary_id) != checked_baries.end())
          continue;
        checked_baries.insert(bary_id);
        el.element = bary_id;
        if (MeshDataMaterialSelector<std::string>::operator()(el) ==
            aggregate_material_id) {
          gel_pockets.push_back(el);
          nb_placed_gel_pockets += 1;
        }
      }
    }
    is_gel_initialized = true;
    std::cout << nb_placed_gel_pockets << " gelpockets placed" << std::endl;
  }

  void initGelPocketPairs() {
    aggregate_material_id = model.getMaterialIndex(aggregate_material);

    Mesh & mesh = this->model.getMesh();
    auto & mesh_facets = mesh.getMeshFacets();
    UInt dim = model.getSpatialDimension();
    //    Element el{_triangle_3, 0, _not_ghost};
    for (auto el_type : model.getMaterial(aggregate_material)
                            .getElementFilter()
                            .elementTypes(dim)) {

      const auto & filter =
          model.getMaterial(aggregate_material).getElementFilter()(el_type);
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

      /// generate the gel pockets
      srand(0.);
      Vector<Real> center(dim);
      std::set<int> checked_baries;
      while (nb_placed_gel_pockets != nb_gel_pockets) {
        /// get a random bary center
        UInt bary_id = rand() % nb_element;
        if (checked_baries.find(bary_id) != checked_baries.end())
          continue;
        checked_baries.insert(bary_id);
        el.element = bary_id;
        if (MeshDataMaterialSelector<std::string>::operator()(el) ==
            aggregate_material_id) {
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
                      connected_element) == aggregate_material_id) {
                gel_pockets.push_back(el);
                gel_pockets.push_back(connected_element);
                nb_placed_gel_pockets += 1;
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
    is_gel_initialized = true;
    std::cout << nb_placed_gel_pockets << " ASR-pocket pairs placed"
              << std::endl;
  }

  UInt operator()(const Element & elem) {

    if (not is_gel_initialized) {
      if (this->gel_pairs) {
        initGelPocketPairs();
      } else {
        initGelPocket();
      }
    }

    UInt temp_index = MeshDataMaterialSelector<std::string>::operator()(elem);
    if (temp_index != aggregate_material_id)
      return temp_index;
    auto iit = gel_pockets.begin();
    auto eit = gel_pockets.end();
    if (std::find(iit, eit, elem) != eit) {
      return model.getMaterialIndex(gel_material);
    }
    return temp_index;
  }

protected:
  SolidMechanicsModel & model;
  std::string gel_material;
  std::vector<Element> gel_pockets;
  UInt nb_gel_pockets;
  UInt nb_placed_gel_pockets;
  std::string aggregate_material{"aggregate"};
  UInt aggregate_material_id{1};
  bool is_gel_initialized{false};
  bool gel_pairs{false};
};

/* --------------------------------------------------------------------------
 */
/* --------------------------------------------------------------------------
 */
using MaterialCohesiveRules = std::map<std::pair<ID, ID>, ID>;

class GelMaterialCohesiveRulesSelector : public MaterialSelector {
public:
  GelMaterialCohesiveRulesSelector(SolidMechanicsModelCohesive & model,
                                   const MaterialCohesiveRules & rules,
                                   std::string gel_material,
                                   const UInt nb_gel_pockets,
                                   std::string aggregate_material = "aggregate",
                                   bool gel_pairs = false)
      : model(model), mesh_data_id("physical_names"), mesh(model.getMesh()),
        mesh_facets(model.getMeshFacets()), dim(model.getSpatialDimension()),
        rules(rules), gel_selector(model, gel_material, nb_gel_pockets,
                                   aggregate_material, gel_pairs),
        default_cohesive(model) {}

  UInt operator()(const Element & element) {
    if (mesh_facets.getSpatialDimension(element.type) == (dim - 1)) {
      const std::vector<Element> & element_to_subelement =
          mesh_facets.getElementToSubelement(element.type, element.ghost_type)(
              element.element);
      // Array<bool> & facets_check = model.getFacetsCheck();

      const Element & el1 = element_to_subelement[0];
      const Element & el2 = element_to_subelement[1];

      ID id1 = model.getMaterial(gel_selector(el1)).getName();
      ID id2 = id1;
      if (el2 != ElementNull) {
        id2 = model.getMaterial(gel_selector(el2)).getName();
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

    return gel_selector(element);
  }

private:
  SolidMechanicsModelCohesive & model;
  ID mesh_data_id;
  const Mesh & mesh;
  const Mesh & mesh_facets;
  UInt dim;
  MaterialCohesiveRules rules;

  GelMaterialSelector gel_selector;
  DefaultMaterialCohesiveSelector default_cohesive;
};

/* ------------------------------------------------------------------------ */
/* Boundary conditions functors */
/* -------------------------------------------------------------------------*/

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
      // if (P < 0)
      //   P = 0.;
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

/* ------------------------------------------------------------------------ */
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
                          const Real ASR_volume_ratio,
                          const std::string group_name,
                          const Real compressibility)
      : model(model), ASR_volume_ratio(ASR_volume_ratio),
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

    // ASR volume computed as AB * thickness (AB * ratio)
    Real ASR_volume = AB_0.norm() * AB_0.norm() * ASR_volume_ratio;
    Real current_volume = AB.norm() * MN.norm();
    current_volume =
        Math::are_float_equal(current_volume, 0.) ? 0 : current_volume;
    Real volume_change = current_volume - ASR_volume;
    Real pressure_change{0};
    if (volume_change < 0) {
      pressure_change = -volume_change / ASR_volume / this->compressibility;
    }
    dual = pressure_change * normal_corrected;
  }

protected:
  SolidMechanicsModel & model;
  const Real ASR_volume_ratio;
  const std::string group_name;
  const Real compressibility;
};
/* ------------------------------------------------------------------- */
class PressureVolumeDependent3D : public BC::Neumann::NeumannFunctor {
public:
  PressureVolumeDependent3D(SolidMechanicsModel & model,
                            const Real ASR_volume_ratio,
                            const std::string group_name,
                            const Real compressibility, const Real multiplier)
      : model(model), ASR_volume_ratio(ASR_volume_ratio),
        group_name(group_name), compressibility(compressibility),
        multiplier(multiplier) {}

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
    if (id < element_ids.size() / 2) {
      normal_corrected *= -1;
      // opposite_facet_nb = element_ids(id + element_ids.size() / 2);
    } else if (id >= element_ids.size()) {
      AKANTU_EXCEPTION("Error in defining side of the cohesive element");
    } else {
      normal_corrected *= 1;
      // opposite_facet_nb = element_ids(id - element_ids.size() / 2);
    }

    // /// compute current area of a gap
    // Vector<UInt> first_facet_nodes = facet_nodes_it[facet_nb];
    // Vector<UInt> second_facet_nodes = facet_nodes_it[opposite_facet_nb];
    // /* corners of a quadrangle consequently
    //         A---M---B
    //         \      |
    //         D--N--C */
    // UInt A, B, C, D;
    // A = first_facet_nodes(0);
    // B = first_facet_nodes(1);
    // C = second_facet_nodes(0);
    // D = second_facet_nodes(1);

    // /// quadrangle's area through diagonals
    // Vector<Real> AC, BD;
    // Vector<Real> A_pos = Vector<Real>(pos_it[A]) +
    // Vector<Real>(disp_it[A]); Vector<Real> B_pos = Vector<Real>(pos_it[B])
    // + Vector<Real>(disp_it[B]); Vector<Real> C_pos =
    // Vector<Real>(pos_it[C]) + Vector<Real>(disp_it[C]); Vector<Real> D_pos
    // = Vector<Real>(pos_it[D]) + Vector<Real>(disp_it[D]); Vector<Real>
    // M_pos = (A_pos + B_pos) * 0.5; Vector<Real> N_pos = (C_pos + D_pos) *
    // 0.5; Vector<Real> MN = M_pos - N_pos; Vector<Real> AB = A_pos - B_pos;
    // Vector<Real> AB_0 = Vector<Real>(pos_it[A]) - Vector<Real>(pos_it[B]);

    // // ASR volume computed as AB * thickness (AB * ratio)
    // Real ASR_volume = AB_0.norm() * AB_0.norm() * ASR_volume_ratio;
    // Real current_volume = AB.norm() * MN.norm();
    // current_volume =
    //     Math::are_float_equal(current_volume, 0.) ? 0 : current_volume;
    // Real volume_change = current_volume - ASR_volume;
    // Real pressure_change{0};
    // if (volume_change < 0) {
    //   pressure_change = -volume_change / ASR_volume /
    //   this->compressibility;
    // }
    // dual = pressure_change * normal_corrected;
    dual = multiplier * normal_corrected;
  }

protected:
  SolidMechanicsModel & model;
  const Real ASR_volume_ratio;
  const std::string group_name;
  const Real compressibility;
  const Real multiplier;
};

/* ------------------------------------------------------------------------ */
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

} // namespace akantu

#endif /* __AKANTU_ASR_TOOLS_HH__ */
