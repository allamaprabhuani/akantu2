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
#include <mesh.hh>
#include <mesh_accessor.hh>
#include <mesh_events.hh>
#include <unordered_set>
/* -------------------------------------------------------------------------- */
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

  template <UInt dim> Real computeSmallestElementSize();

  // /// apply homogeneous temperature on Solid Mechanics Model
  // void applyTemperatureFieldToSolidmechanicsModel(const Real & temperature);

  /// compute increase in gel strain within 1 timestep
  Real computeDeltaGelStrainThermal(const Real delta_time, const Real k,
                                    const Real activ_energy, const Real R,
                                    const Real T,
                                    Real & amount_reactive_particles,
                                    const Real saturation_const);

  /// compute linear increase in gel strain
  Real computeDeltaGelStrainLinear(const Real delta_time, const Real k);

  /// insert single cohesive element by the coordinate of its center
  void insertCohElemByCoord(const Vector<Real> & position);

  /// insert multiple cohesive elements by the limiting box
  void insertCohElemByLimits(const Matrix<Real> & insertion_limits,
                             std::string coh_mat_name);

  /// insert multiple cohesive elements by the limiting box
  void insertCohElemRandomly(const UInt & nb_coh_elem, std::string coh_mat_name,
                             std::string matrix_mat_name);

  /// insert up to 3 facets pair based on the coord of the central one
  const Array<Element>
  insertCohElOrFacetsByCoord(const Vector<Real> & position,
                             bool add_neighbors = true,
                             bool only_double_facets = false);

  /// on elements added for asr-tools
  void onElementsAdded(const Array<Element> & elements,
                       const NewElementsEvent & element_event);

  void onNodesAdded(const Array<UInt> & new_nodes,
                    const NewNodesEvent & nodes_event);

  /// apply self-weight force
  void applyBodyForce();

  /// apply delta u on nodes
  void applyDeltaU(Real delta_u);

  /// apply eigenstrain on gel material
  void applyGelStrain(const Matrix<Real> & prestrain);
  /* ------------------------------------------------------------------------ */
  /// RVE part

  /// apply boundary contions based on macroscopic deformation gradient
  virtual void
  applyBoundaryConditionsRve(const Matrix<Real> & displacement_gradient);

  /// apply homogeneous temperature field from the macroscale level to the RVEs
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

  /// dump the RVE
  void dumpRve();

  /// compute average stress in the RVE
  void homogenizeStressField(Matrix<Real> & stress);
protected:
    /// storing nodal fields before tests
  void storeNodalFields();

  /// restoring nodal fields before tests
  void restoreNodalFields();

  /// restoring internal fields after tests
  void restoreInternalFields();


private:
  /// find the corner nodes
  void findCornerNodes();

  /// perform virtual testing
  void performVirtualTesting(const Matrix<Real> & H,
                             Matrix<Real> & eff_stresses,
                             Matrix<Real> & eff_strains, const UInt test_no);

  /* ------------------------------------------------------------------------
   */
public:
  // Accessors
  inline Array<std::tuple<UInt, UInt>> getNodePairs() const {
    return node_pairs;
  }

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

  /// lower limit for stresses
  Array<Real> stress_limit;

  /// dump counter
  UInt nb_dumps;

  /// booleans for applying delta u
  bool doubled_facets_ready;
  bool doubled_nodes_ready;

  // array of tuples to store nodes pairs:1st- is the one on the upper facet
  Array<std::tuple<UInt, UInt>> node_pairs;

  // arrays to store nodal values during virtual tests
  Array<Real> disp_stored;
  Array<Real> ext_force_stored;
  Array<bool> boun_stored;
};

/* -------------------------------------------------------------------------- */
/* ASR material selector                                                      */
/* -------------------------------------------------------------------------- */
class GelMaterialSelector : public MeshDataMaterialSelector<std::string> {
public:
  GelMaterialSelector(SolidMechanicsModel & model, const Real box_size,
                      const std::string & gel_material,
                      const UInt nb_gel_pockets,
                      std::string aggregate_material = "aggregate",
                      Real /*tolerance*/ = 0.)
      : MeshDataMaterialSelector<std::string>("physical_names", model),
        model(model), gel_material(gel_material),
        nb_gel_pockets(nb_gel_pockets), nb_placed_gel_pockets(0),
        box_size(box_size), aggregate_material(aggregate_material) {}

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
      UInt placed_gel_pockets = 0;
      std::set<int> checked_baries;
      while (placed_gel_pockets != nb_gel_pockets) {
        /// get a random bary center
        UInt bary_id = rand() % nb_element;
        if (checked_baries.find(bary_id) != checked_baries.end())
          continue;
        checked_baries.insert(bary_id);
        el.element = bary_id;
        if (MeshDataMaterialSelector<std::string>::operator()(el) ==
            aggregate_material_id) {
          gel_pockets.push_back(el);
          placed_gel_pockets += 1;
        }
      }
    }
    is_gel_initialized = true;
  }

  UInt operator()(const Element & elem) {
    /// variables for parallel execution
    auto && comm = akantu::Communicator::getWorldCommunicator();
    auto prank = comm.whoAmI();

    if (not is_gel_initialized)
      initGelPocket();

    UInt temp_index = MeshDataMaterialSelector<std::string>::operator()(elem);
    if (temp_index != aggregate_material_id)
      return temp_index;
    auto iit = gel_pockets.begin();
    auto eit = gel_pockets.end();
    if (std::find(iit, eit, elem) != eit) {
      nb_placed_gel_pockets += 1;
      if (prank == 0)
        std::cout << nb_placed_gel_pockets << " gelpockets placed" << std::endl;
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
  Real box_size;
  std::string aggregate_material{"aggregate"};
  UInt aggregate_material_id{1};
  bool is_gel_initialized{false};
};

/* -------------------------------------------------------------------------- */
/* Boundary conditions functors */
/* -------------------------------------------------------------------------- */

class Pressure : public BC::Neumann::NeumannFunctor {
public:
  Pressure(SolidMechanicsModel & model,
           const Array<akantu::Real> & pressure_on_qpoint,
           const Array<akantu::Real> & quad_coords)
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

/* -------------------------------------------------------------------------- */

class DeltaU : public BC::Dirichlet::DirichletFunctor {
public:
  DeltaU(const SolidMechanicsModel & model, const Real delta_u,
         const Array<std::tuple<UInt, UInt>> & node_pairs)
      : model(model), delta_u(delta_u), node_pairs(node_pairs) {
    displacement = model.getDisplacement();
  }

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
