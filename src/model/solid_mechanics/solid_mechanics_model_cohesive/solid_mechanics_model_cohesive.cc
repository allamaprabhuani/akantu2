/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model_cohesive.hh"
#include "aka_iterators.hh"
#include "cohesive_element_inserter.hh"
#include "element_synchronizer.hh"
#include "facet_synchronizer.hh"
#include "fe_engine_template.hh"
#include "global_ids_updater.hh"
#include "integrator_gauss.hh"
#include "material_cohesive.hh"
#include "mesh_accessor.hh"
#include "mesh_global_data_updater.hh"
#include "parser.hh"
#include "shape_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

namespace akantu {

class CohesiveMeshGlobalDataUpdater : public MeshGlobalDataUpdater {
public:
  CohesiveMeshGlobalDataUpdater(SolidMechanicsModelCohesive & model)
      : model(model), mesh(model.getMesh()),
        global_ids_updater(model.getMesh(), model.cohesive_synchronizer.get()) {
  }

  /* ------------------------------------------------------------------------ */
  std::tuple<UInt, UInt>
  updateData(NewNodesEvent & nodes_event,
             NewElementsEvent & elements_event) override {
    auto * cohesive_nodes_event =
        dynamic_cast<CohesiveNewNodesEvent *>(&nodes_event);
    if (cohesive_nodes_event == nullptr) {
      return std::make_tuple(nodes_event.getList().size(),
                             elements_event.getList().size());
    }

    /// update nodes type
    auto & new_nodes = cohesive_nodes_event->getList();
    auto & old_nodes = cohesive_nodes_event->getOldNodesList();

    auto local_nb_new_nodes = new_nodes.size();
    auto nb_new_nodes = local_nb_new_nodes;

    if (mesh.isDistributed()) {
      MeshAccessor mesh_accessor(mesh);
      auto & nodes_flags = mesh_accessor.getNodesFlags();
      auto nb_old_nodes = nodes_flags.size();
      nodes_flags.resize(nb_old_nodes + local_nb_new_nodes);

      for (auto && [old_node, new_node] : zip(old_nodes, new_nodes)) {
        nodes_flags(new_node) = nodes_flags(old_node);
      }

      model.updateCohesiveSynchronizers(elements_event);
      nb_new_nodes = global_ids_updater.updateGlobalIDs(new_nodes.size());
    }

    auto nb_new_elements = elements_event.getList().size();
    const auto & comm = mesh.getCommunicator();
    comm.allReduce(nb_new_elements, SynchronizerOperation::_sum);

    if (nb_new_elements > 0) {
      mesh.sendEvent(elements_event);
    }

    if (nb_new_nodes > 0) {
      mesh.sendEvent(nodes_event);
    }

    return std::make_tuple(nb_new_nodes, nb_new_elements);
  }

private:
  SolidMechanicsModelCohesive & model;
  Mesh & mesh;
  GlobalIdsUpdater global_ids_updater;
};

/* -------------------------------------------------------------------------- */
SolidMechanicsModelCohesive::SolidMechanicsModelCohesive(
    Mesh & mesh, Int dim, const ID & id,
    const std::shared_ptr<DOFManager> & dof_manager)
    : SolidMechanicsModel(mesh, dim, id, dof_manager,
                          ModelType::_solid_mechanics_model_cohesive),
      tangents("tangents", id), facet_stress("facet_stress", id),
      facet_material("facet_material", id) {
  AKANTU_DEBUG_IN();

  registerFEEngineObject<MyFEEngineCohesiveType>("CohesiveFEEngine", mesh,
                                                 Model::spatial_dimension);

  auto && tmp_material_selector =
      std::make_shared<DefaultMaterialCohesiveSelector>(*this);

  tmp_material_selector->setFallback(this->getConstitutiveLawSelector());
  this->setConstitutiveLawSelector(tmp_material_selector);

  this->mesh.registerDumper<DumperParaview>("cohesive elements", id);
  this->mesh.addDumpMeshToDumper("cohesive elements", mesh,
                                 Model::spatial_dimension, _not_ghost,
                                 _ek_cohesive);

  if (this->mesh.isDistributed()) {
    /// create the distributed synchronizer for cohesive elements
    this->cohesive_synchronizer = std::make_unique<ElementSynchronizer>(
        mesh, "cohesive_distributed_synchronizer");

    auto & synchronizer = mesh.getElementSynchronizer();
    this->cohesive_synchronizer->split(synchronizer, [](auto && el) {
      return Mesh::getKind(el.type) == _ek_cohesive;
    });

    this->registerSynchronizer(*cohesive_synchronizer,
                               SynchronizationTag::_constitutive_law_id);
    this->registerSynchronizer(*cohesive_synchronizer,
                               SynchronizationTag::_smm_stress);
    this->registerSynchronizer(*cohesive_synchronizer,
                               SynchronizationTag::_smm_boundary);
  }

  this->inserter = std::make_unique<CohesiveElementInserter>(
      this->mesh, id + ":cohesive_element_inserter");

  registerFEEngineObject<MyFEEngineFacetType>(
      "FacetsFEEngine", mesh.getMeshFacets(), Model::spatial_dimension - 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::setTimeStep(Real time_step,
                                              const ID & solver_id) {
  SolidMechanicsModel::setTimeStep(time_step, solver_id);
  this->mesh.getDumper("cohesive elements").setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initFullImpl(const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  const auto & smmc_options =
      aka::as_type<SolidMechanicsModelCohesiveOptions>(options);

  this->is_extrinsic = smmc_options.is_extrinsic;

  inserter->setIsExtrinsic(is_extrinsic);

  if (mesh.isDistributed()) {
    auto & mesh_facets = inserter->getMeshFacets();
    auto & synchronizer =
        aka::as_type<FacetSynchronizer>(mesh_facets.getElementSynchronizer());

    // synchronizeGhostFacetsConnectivity();

    /// create the facet synchronizer for extrinsic simulations
    if (is_extrinsic) {
      facet_stress_synchronizer = std::make_unique<ElementSynchronizer>(
          synchronizer, id + ":facet_stress_synchronizer");
      facet_stress_synchronizer->swapSendRecv();
      this->registerSynchronizer(*facet_stress_synchronizer,
                                 SynchronizationTag::_smmc_facets_stress);
    }
  }

  MeshAccessor mesh_accessor(mesh);
  mesh_accessor.registerGlobalDataUpdater(
      std::make_unique<CohesiveMeshGlobalDataUpdater>(*this));

  auto && [section, is_empty] = this->getParserSection();

  if (not is_empty) {
    auto inserter_section =
        section.getSubSections(ParserType::_cohesive_inserter);
    if (inserter_section.begin() != inserter_section.end()) {
      inserter->parseSection(*inserter_section.begin());
    }
  }

  SolidMechanicsModel::initFullImpl(options);

  AKANTU_DEBUG_OUT();
} // namespace akantu

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initConstitutiveLaws() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  instantiateMaterials();

  auto & material_selector = this->getConstitutiveLawSelector();

  // set the facet information in the material in case of dynamic insertion
  // to know what material to call for stress checks
  const Mesh & mesh_facets = inserter->getMeshFacets();
  facet_material.initialize(
      mesh_facets, _spatial_dimension = spatial_dimension - 1,
      _with_nb_element = true,
      _default_value =
          DefaultMaterialCohesiveSelector::getDefaultCohesiveMaterial(*this));

  for_each_element(
      mesh_facets,
      [&](auto && element) {
        auto mat_index = material_selector(element);
        if (not mat_index) {
          return;
        }
        auto & mat =
            aka::as_type<MaterialCohesive>(this->getConstitutiveLaw(mat_index));

        facet_material(element) = mat_index;
        if (is_extrinsic) {
          mat.addFacet(element);
        }
      },
      _spatial_dimension = spatial_dimension - 1, _ghost_type = _not_ghost);

  SolidMechanicsModel::initConstitutiveLaws();

  auto & initial_nodes = mesh.getNodalData<Idx>("initial_nodes_match");
  initial_nodes.resize(mesh.getNbNodes());
  for (auto && [node, init_node] : enumerate(initial_nodes)) {
    init_node = node;
  }

  lambda = std::make_unique<Array<Real>>(0, spatial_dimension, "lambda");
  lambda_blocked_dofs = std::make_unique<Array<bool>>(0, spatial_dimension, "lambda_blocked_dofs");

  if (lambda) {
    mesh.getElementalData<Idx>("initial_nodes_connectivities");
    mesh.getElementalData<Idx>("lambda_connectivities");
    auto & nodes_to_lambda = mesh.getNodalData<Idx>("nodes_to_lambda");
    nodes_to_lambda.resize(mesh.getNbNodes(), -1);
  }

  if (is_extrinsic) {
    this->initAutomaticInsertion();
  } else {
    this->insertIntrinsicElements();
  }

  if (lambda) {
    auto & dof_manager = this->getDOFManager();
    if (!dof_manager.hasDOFs("lambda")) {
        dof_manager.registerDOFs("lambda", *this->lambda, _dst_generic);
        dof_manager.registerBlockedDOFs("lambda", *this->lambda_blocked_dofs);
    }
  }

  AKANTU_DEBUG_OUT();
} // namespace akantu

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model,basically it  pre-compute the shapes, shapes derivatives
 * and jacobian
 */
void SolidMechanicsModelCohesive::initModel() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::initModel();

  /// add cohesive type connectivity
  ElementType type = _not_defined;
  for (auto && type_ghost : ghost_types) {
    for (const auto & tmp_type :
         mesh.elementTypes(spatial_dimension, type_ghost)) {
      const auto & connectivity = mesh.getConnectivity(tmp_type, type_ghost);
      if (connectivity.empty()) {
        continue;
      }

      type = tmp_type;
      auto type_facet = Mesh::getFacetType(type);
      auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);
      AKANTU_DEBUG_ASSERT(Mesh::getKind(type_cohesive) == _ek_cohesive,
                          "The element type " << type_cohesive
                                              << " is not a cohesive type");
      mesh.addConnectivityType(type_cohesive, type_ghost);
    }
  }
  AKANTU_DEBUG_ASSERT(type != _not_defined, "No elements in the mesh");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::insertIntrinsicElements() {
  AKANTU_DEBUG_IN();
  inserter->insertIntrinsicElements();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  this->inserter->limitCheckFacets();
  this->updateFacetStressSynchronizer();
  this->resizeFacetStress();

  /// compute normals on facets
  this->computeNormals();

  this->initStressInterpolation();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::updateAutomaticInsertion() {
  AKANTU_DEBUG_IN();

  this->inserter->limitCheckFacets();
  this->updateFacetStressSynchronizer();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::initStressInterpolation() {
  Mesh & mesh_facets = inserter->getMeshFacets();

  /// compute quadrature points coordinates on facets
  Array<Real> & position = mesh.getNodes();

  ElementTypeMapArray<Real> quad_facets("quad_facets", id);
  quad_facets.initialize(mesh_facets, _nb_component = Model::spatial_dimension,
                         _spatial_dimension = Model::spatial_dimension - 1);

  getFEEngine("FacetsFEEngine")
      .interpolateOnIntegrationPoints(position, quad_facets);

  /// compute elements quadrature point positions and build
  /// element-facet quadrature points data structure
  ElementTypeMapArray<Real> elements_quad_facets("elements_quad_facets", id);

  elements_quad_facets.initialize(
      mesh, _nb_component = Model::spatial_dimension,
      _spatial_dimension = Model::spatial_dimension);

  for (auto elem_gt : ghost_types) {
    for (const auto & type :
         mesh.elementTypes(Model::spatial_dimension, elem_gt)) {
      auto nb_element = mesh.getNbElement(type, elem_gt);
      if (nb_element == 0) {
        continue;
      }

      /// compute elements' quadrature points and list of facet
      /// quadrature points positions by element
      const auto & facet_to_element =
          mesh_facets.getSubelementToElement(type, elem_gt);
      auto & el_q_facet = elements_quad_facets(type, elem_gt);

      auto facet_type = Mesh::getFacetType(type);
      auto nb_quad_per_facet =
          getFEEngine("FacetsFEEngine").getNbIntegrationPoints(facet_type);
      auto nb_facet_per_elem = facet_to_element.getNbComponent();

      // small hack in the loop to skip boundary elements, they are silently
      // initialized to NaN to see if this causes problems
      el_q_facet.resize(nb_element * nb_facet_per_elem * nb_quad_per_facet,
                        std::numeric_limits<Real>::quiet_NaN());

      for (auto && data :
           zip(make_view(facet_to_element),
               make_view(el_q_facet, spatial_dimension, nb_quad_per_facet))) {
        const auto & global_facet = std::get<0>(data);
        auto & el_q = std::get<1>(data);

        if (global_facet == ElementNull) {
          continue;
        }

        auto && quad_f =
            make_view(quad_facets(global_facet.type, global_facet.ghost_type),
                      spatial_dimension, nb_quad_per_facet)
                .begin()[global_facet.element];

        el_q = quad_f;
      }
    }
  }

  /// loop over non cohesive materials
  this->for_each_constitutive_law([&](auto && material) {
    if (aka::is_of_type<MaterialCohesive>(material)) {
      return;
    }
    /// initialize the interpolation function
    material.initElementalFieldInterpolation(elements_quad_facets);
  });

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  // f_int += f_int_cohe
  this->for_each_constitutive_law([&](auto && material) {
    try {
      auto & mat = aka::as_type<MaterialCohesive>(material);
      mat.computeTraction(_not_ghost);
    } catch (std::bad_cast & bce) {
    }
  });

//  ArrayPrintHelper<true>::ArrayPrintHelper::print_content<bool>(*(this->blocked_dofs),std::cout,0);
  SolidMechanicsModel::assembleInternalForces();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::computeNormals() {
  AKANTU_DEBUG_IN();

  Mesh & mesh_facets = this->inserter->getMeshFacets();
  this->getFEEngine("FacetsFEEngine")
      .computeNormalsOnIntegrationPoints(_not_ghost);

  /**
   *  @todo store tangents while computing normals instead of
   *  recomputing them as follows:
   */
  /* ------------------------------------------------------------------------ */
  UInt tangent_components =
      Model::spatial_dimension * (Model::spatial_dimension - 1);

  tangents.initialize(mesh_facets, _nb_component = tangent_components,
                      _spatial_dimension = Model::spatial_dimension - 1);

  for (auto facet_type :
       mesh_facets.elementTypes(Model::spatial_dimension - 1)) {
    const Array<Real> & normals =
        this->getFEEngine("FacetsFEEngine")
            .getNormalsOnIntegrationPoints(facet_type);

    Array<Real> & tangents = this->tangents(facet_type);

    Math::compute_tangents(normals, tangents);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::interpolateStress() {
  ElementTypeMapArray<Real> by_elem_result("temporary_stress_by_facets", id);

  this->for_each_constitutive_law([&](auto && material) {
    if (not aka::is_of_type<MaterialCohesive>(material)) {
      /// interpolate stress on facet quadrature points positions
      material.interpolateStressOnFacets(facet_stress, by_elem_result);
    }
  });

  this->synchronize(SynchronizationTag::_smmc_facets_stress);
}

/* -------------------------------------------------------------------------- */
UInt SolidMechanicsModelCohesive::checkCohesiveStress() {
  AKANTU_DEBUG_IN();

  if (not is_extrinsic) {
    AKANTU_EXCEPTION(
        "This function can only be used for extrinsic cohesive elements");
  }

  interpolateStress();

  this->for_each_constitutive_law([&](auto && material) {
    if (aka::is_of_type<MaterialCohesive>(material)) {
      /// check which not ghost cohesive elements are to be created
      auto & mat_cohesive = aka::as_type<MaterialCohesive>(material);
      mat_cohesive.checkInsertion();
    }
  });

  /// insert cohesive elements
  auto nb_new_elements = inserter->insertElements();

  AKANTU_DEBUG_OUT();

  return nb_new_elements;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onElementsAdded(
    const Array<Element> & element_list, const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::onElementsAdded(element_list, event);

  if (lambda) {
    auto & initial_nodes_connectivities =
        mesh.getElementalData<Idx>("initial_nodes_connectivities");
    auto & lambda_connectivities =
        mesh.getElementalData<Idx>("lambda_connectivities");

    for (auto ghost_type : ghost_types) {
      for (auto type : mesh.elementTypes(_element_kind = _ek_cohesive,
                       _ghost_type = ghost_type)) {
        auto size = mesh.getConnectivity(type, ghost_type).size();
        if (not initial_nodes_connectivities.exists(type, ghost_type)) {
          auto underlying_type = Mesh::getFacetType(type);
          initial_nodes_connectivities.alloc(
              size, Mesh::getNbNodesPerElement(underlying_type), type,
              ghost_type, -1);
          lambda_connectivities.alloc(
              size, Mesh::getNbNodesPerElement(underlying_type), type,
              ghost_type, -1);

        } else {
          initial_nodes_connectivities(type, ghost_type).resize(size, -1);
          lambda_connectivities(type, ghost_type).resize(size, -1);
        }
      }
    }

    // Set some connectivities, it will be corrected at on nodes added
    for (const auto & el : element_list) {
      if (Mesh::getKind(el.type) != _ek_cohesive) {
        continue;
      }

      auto && conn = mesh.getConnectivity(el);
      auto && iconn = initial_nodes_connectivities.get(el);
      iconn = conn.block(0, 0, iconn.size(), 1);
    }
  }

  if (is_extrinsic) {
    resizeFacetStress();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onNodesAdded(const Array<Idx> & new_nodes,
                                               const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::onNodesAdded(new_nodes, event);

  auto & initial_nodes = mesh.getNodalData<Idx>("initial_nodes_match");
  auto old_max_nodes = initial_nodes.size();
  initial_nodes.resize(mesh.getNbNodes());

  const auto * cohesive_event =
      dynamic_cast<const CohesiveNewNodesEvent *>(&event);
  if (cohesive_event == nullptr) {
    for (auto && node : new_nodes) {
      initial_nodes(node) = node;
    }
    return;
  }

  auto old_nodes = cohesive_event->getOldNodesList();

  // getting a corrected old_nodes for multiple doubling
  for (auto & onode : old_nodes) {
    while (onode >= old_max_nodes) {
      auto nnode = new_nodes.find(onode);
      AKANTU_DEBUG_ASSERT(nnode != -1,
                          "If the old node is bigger than old_max_nodes it "
                          "should also be a new node");
      onode = nnode;
    }
  }

  auto copy = [&new_nodes, &old_nodes](auto & arr) {
    auto it = make_view(arr, arr.getNbComponent()).begin();
    for (auto [new_node, old_node] : zip(new_nodes, old_nodes)) {
      it[new_node] = it[old_node];
    }
  };

  copy(*displacement);
  copy(*blocked_dofs);

  if (velocity) {
    copy(*velocity);
  }

  if (acceleration) {
    copy(*acceleration);
  }

  if (current_position) {
    copy(*current_position);
  }

  if (previous_displacement) {
    copy(*previous_displacement);
  }

  if (displacement_increment) {
    copy(*displacement_increment);
  }

  copy(getDOFManager().getSolution("displacement"));

  copy(initial_nodes);  this->allocNodalField(this->blocked_dofs, spatial_dimension, "blocked_dofs");


  // correct connectivities
  if (lambda) {
    auto & initial_nodes_connectivities =
        mesh.getElementalData<Idx>("initial_nodes_connectivities");
    auto & lambda_connectivities =
        mesh.getElementalData<Idx>("lambda_connectivities");

    auto & nodes_to_lambda = mesh.getNodalData<Idx>("nodes_to_lambda");

    auto lambda_id = lambda->size();

    for (auto ghost_type : ghost_types) {
      for (auto type : initial_nodes_connectivities.elementTypes(
               _element_kind = _ek_cohesive, _ghost_type = ghost_type)) {

        auto & initial_nodes_connectivity =
            initial_nodes_connectivities(type, ghost_type);
        auto & lambda_connectivity = lambda_connectivities(type, ghost_type);

        for (auto && [conn, lambda_conn] :
             zip(make_view(initial_nodes_connectivity),
                 make_view(lambda_connectivity))) {
          conn = initial_nodes(conn);
          auto & ntl = nodes_to_lambda(conn);
          if (ntl == -1) {
            ntl = lambda_id;
            ++lambda_id;
          }

          lambda_conn = ntl;
        }
      }
    }

    lambda->resize(lambda_id, 0.);
    lambda_blocked_dofs->resize(lambda_id, false);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::afterSolveStep(bool converged) {
  AKANTU_DEBUG_IN();

  /*
   * This is required because the Cauchy stress is the stress measure that
   * is used to check the insertion of cohesive elements
   */
  if (converged) {
    this->for_each_constitutive_law([](auto && mat) {
      if (mat.isFiniteDeformation()) {
        mat.computeAllCauchyStresses(_not_ghost);
      }
    });
  }

  SolidMechanicsModel::afterSolveStep(converged);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions SolidMechanicsModelCohesive::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options = SolidMechanicsModel::getDefaultSolverOptions(type);

  if(lambda)
  {
      switch (type) {
      case TimeStepSolverType::_dynamic_lumped: {
          options.non_linear_solver_type = NonLinearSolverType::_lumped;
          options.integration_scheme_type["lambda"] =
                  IntegrationSchemeType::_central_difference;
          options.solution_type["lambda"] = IntegrationScheme::_acceleration;
          break;
      }
      case TimeStepSolverType::_static: {
          options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
          options.integration_scheme_type["lambda"] =
                  IntegrationSchemeType::_pseudo_time;
          options.solution_type["lambda"] = IntegrationScheme::_not_defined;
          break;
      }
      case TimeStepSolverType::_dynamic: {
          if (this->method == _explicit_consistent_mass) {
              options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
              options.integration_scheme_type["lambda"] =
                      IntegrationSchemeType::_central_difference;
              options.solution_type["lambda"] = IntegrationScheme::_acceleration;
          } else {
              options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
              options.integration_scheme_type["lambda"] =
                      IntegrationSchemeType::_trapezoidal_rule_2;
              options.solution_type["lambda"] = IntegrationScheme::_displacement;
          }
          break;
      }
      default:
          AKANTU_EXCEPTION(type << " is not a valid time step solver type");
      }

  }
  return options;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::printself(std::ostream & stream,
                                            int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "SolidMechanicsModelCohesive ["
         << "\n";
  SolidMechanicsModel::printself(stream, indent + 2);
  stream << space << "]\n";
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::resizeFacetStress() {
  AKANTU_DEBUG_IN();

  this->facet_stress.initialize(getFEEngine("FacetsFEEngine"),
                                _nb_component =
                                    2 * spatial_dimension * spatial_dimension,
                                _spatial_dimension = spatial_dimension - 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::addDumpGroupFieldToDumper(
    const std::string & dumper_name, const std::string & field_id,
    const std::string & group_name, ElementKind element_kind,
    bool padding_flag) {
  AKANTU_DEBUG_IN();

  Int spatial_dimension = Model::spatial_dimension;
  ElementKind _element_kind = element_kind;
  if (dumper_name == "cohesive elements") {
    _element_kind = _ek_cohesive;
  } else if (dumper_name == "facets") {
    spatial_dimension = Model::spatial_dimension - 1;
  }

  SolidMechanicsModel::addDumpGroupFieldToDumper(dumper_name, field_id,
                                                 group_name, spatial_dimension,
                                                 _element_kind, padding_flag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelCohesive::onDump() {
  this->flattenAllRegisteredInternals(_ek_cohesive);
  SolidMechanicsModel::onDump();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
