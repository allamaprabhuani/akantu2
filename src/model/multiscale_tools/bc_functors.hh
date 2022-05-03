/**
 * @file   bc_functors.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue Jan 16 10:26:53 2014
 * @update Thursday Mar 24  2022
 * @brief  functors to apply internal loading (e.g. pressure in cracks)
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
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

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
