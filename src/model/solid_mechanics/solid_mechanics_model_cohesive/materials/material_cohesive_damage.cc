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
#include "material_cohesive_damage.hh"
#include "dof_synchronizer.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "fe_engine_template.hh"
#include "integrator_gauss.hh"
#include "shape_cohesive.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialCohesiveDamage<dim>::MaterialCohesiveDamage(SolidMechanicsModel & model,
                                                    const ID & id)
    : MaterialCohesive(model, id),
      czm_damage(registerInternal<Real, CohesiveInternalField>("czm_damage", 1)),
      lambda(registerInternal<Real, CohesiveInternalField>("lambda", dim)){
  AKANTU_DEBUG_IN();

  this->registerParam("k", k, Real(0.), _pat_parsable | _pat_readable,
                      "Beta parameter");

  this->registerParam("G_c", G_c, Real(0.), _pat_parsable | _pat_readable,
                      "Mode I fracture energy");


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialCohesiveDamage<dim>::initMaterial() {
    AKANTU_DEBUG_IN();

    MaterialCohesive::initMaterial();
    //  lambda.initialize(dim);

    //  const auto & mesh_facets = model->getMeshFacets();
    //  for (const auto & type_facet : mesh_facets.elementTypes(dim - 1)) {
    //    auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);

    //    const auto & facet_filter_array = facet_filter(type_facet);
    //    const auto & lambda_array = lambda(type_cohesive);

    //    auto nb_quad_facet =
    //        model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);

    ////    for (auto && [facet, lda] :
    ////         zip(facet_filter_array, lambda_array)) {

    //////        for (Int q = 0; q < nb_quad_facet; ++q) {
    //////            auto current_quad = facet * nb_quad_facet + q;
    //////        }
    ////        std::cout << " lda = " << lda << std::endl;
    ////    }

    ////    auto & dof_manager = this->model->getDOFManager();

    ////    std::unique_ptr<Array<Real>> toto;

    ////    dof_manager.registerDOFs("lambda", *toto, _dst_generic);

    //  }

//    GhostType ghost_type = _not_ghost;

//    for (auto type : getElementFilter().elementTypes(spatial_dimension,
//                                                     ghost_type, _ek_cohesive)) {
//        auto & elem_filter = getElementFilter(type, ghost_type);
//        auto nb_element = elem_filter.size();
//        if (nb_element == 0) {
//            continue;
//        }

//        auto & lambda = lambdas(type, ghost_type);

//        auto lambda_it = lambda.begin(dim, 1);

//        auto nb_quadrature_points =
//                fem_cohesive.getNbIntegrationPoints(type, ghost_type);

//        for (Int el = 0; el < nb_element; ++el) {
//            auto current_quad = elem_filter(el) * nb_quadrature_points;

//            for (Int q = 0; q < nb_quadrature_points; ++q, ++lambda_it) {

//                //        std::cout << "lda = ",*lambda_it << std::endl;
//            }
//        }

//        auto & dof_manager = this->model->getDOFManager();
//        dof_manager.registerDOFs("lambdas", lambdas, _dst_generic);


//    }
    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

//  auto & internal_force = const_cast<Array<Real> &>(model->getInternalForce());

//  for (auto type : getElementFilter().elementTypes(spatial_dimension,
//                                                   ghost_type, _ek_cohesive)) {
//    auto & elem_filter = getElementFilter(type, ghost_type);
//    auto nb_element = elem_filter.size();
//    if (nb_element == 0) {
//      continue;
//    }

//    const auto & shapes = fem_cohesive.getShapes(type, ghost_type);

//    auto size_of_shapes = shapes.getNbComponent();
//    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
//    auto nb_quadrature_points =
//        fem_cohesive.getNbIntegrationPoints(type, ghost_type);

//  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::computeLambdaOnQuad(const Array<Real> & lambda_nodes,
                                      Array<Real> & lambda_quad, ElementType type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem_cohesive =
      this->model->getFEEngineClass<MyFEEngineCohesiveType>("CohesiveFEEngine");

  //// TODO : THIS IS WRONG, SHOULD BE CORRECTED
  tuple_dispatch<ElementTypes_t<_ek_cohesive>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        fem_cohesive.getShapeFunctions()
            .interpolateOnIntegrationPoints<type,
                                            CohesiveReduceFunctionOpening>(
                lambda_nodes, lambda_quad, spatial_dimension, ghost_type,
                this->getElementFilter(type, ghost_type));
      },
      type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::computeTraction(ElementType el_type,
                                                  GhostType ghost_type) {

    computeLambdaOnQuad(model->getLambda(),lambda(el_type, ghost_type),el_type,ghost_type);
    for (auto && args : getArguments(el_type, ghost_type)) {
      this->computeTractionOnQuad(args);
    }
}

/* -------------------------------------------------------------------------- */
template class MaterialCohesiveDamage<1>;
template class MaterialCohesiveDamage<2>;
template class MaterialCohesiveDamage<3>;
const bool material_is_alocated_cohesive_damage [[maybe_unused]] =
    instantiateMaterial<MaterialCohesiveDamage>("cohesive_damage");

} // namespace akantu
