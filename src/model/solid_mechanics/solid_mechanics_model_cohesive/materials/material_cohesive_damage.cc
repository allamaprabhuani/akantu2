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
#include "fe_engine_template.hh"
#include "integrator_gauss.hh"
#include "shape_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
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
      lambda(registerInternal<Real, CohesiveInternalField>("lambda",
                                                           spatial_dimension)),
      err_openings(registerInternal<Real, CohesiveInternalField>(
          "err_openings", spatial_dimension)),
      czm_damage(
          registerInternal<Real, CohesiveInternalField>("czm_damage", 1)) {
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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & internal_force = const_cast<Array<Real> &>(model->getInternalForce());

  for (auto type : getElementFilter().elementTypes(spatial_dimension,
                                                   ghost_type, _ek_cohesive)) {
    auto & elem_filter = getElementFilter(type, ghost_type);
    auto nb_element = elem_filter.size();
    if (nb_element == 0) {
      continue;
    }

    const auto & shapes = fem_cohesive.getShapes(type, ghost_type);
    auto & traction = tractions(type, ghost_type);
    auto & err_opening = err_openings(type, ghost_type);

    auto size_of_shapes = shapes.getNbComponent();
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points =
        fem_cohesive.getNbIntegrationPoints(type, ghost_type);

    /// compute @f$t_i N_a@f$

    auto traction_cpy = std::make_shared<Array<Real>>(
        nb_element * nb_quadrature_points, spatial_dimension * size_of_shapes);

    auto err_opening_cpy = std::make_shared<Array<Real>>(
        nb_element * nb_quadrature_points, spatial_dimension * size_of_shapes);

    auto traction_it = traction.begin(spatial_dimension, 1);
    auto err_opening_it = err_opening.begin(spatial_dimension, 1);

    auto shapes_filtered_begin = shapes.begin(1, size_of_shapes);
    auto traction_cpy_it =
        traction_cpy->begin(spatial_dimension, size_of_shapes);
    auto err_opening_cpy_it =
        err_opening_cpy->begin(spatial_dimension, size_of_shapes);

    for (Int el = 0; el < nb_element; ++el) {
      auto current_quad = elem_filter(el) * nb_quadrature_points;

      for (Int q = 0; q < nb_quadrature_points; ++q, ++traction_it,
               ++err_opening_it, ++current_quad, ++traction_cpy_it,
               ++err_opening_cpy_it) {

        auto && shapes_filtered = shapes_filtered_begin[current_quad];

        *traction_cpy_it = (*traction_it) * shapes_filtered;
        *err_opening_cpy_it = (*err_opening_it) * shapes_filtered;
      }
    }

    /**
     * compute @f$\int t \cdot N\, dS@f$ by  @f$ \sum_q \mathbf{N}^t
     * \mathbf{t}_q \overline w_q J_q@f$
     */
    auto partial_int_t_N = std::make_shared<Array<Real>>(
        nb_element, spatial_dimension * size_of_shapes, "int_t_N");

    fem_cohesive.integrate(*traction_cpy, *partial_int_t_N,
                           spatial_dimension * size_of_shapes, type, ghost_type,
                           elem_filter);

    auto int_t_N = std::make_shared<Array<Real>>(
        nb_element, 2 * spatial_dimension * size_of_shapes, "int_t_N");

    auto * int_t_N_val = int_t_N->data();
    auto * partial_int_t_N_val = partial_int_t_N->data();
    for (Int el = 0; el < nb_element; ++el) {
      std::copy_n(partial_int_t_N_val, size_of_shapes * spatial_dimension,
                  int_t_N_val);
      std::copy_n(partial_int_t_N_val, size_of_shapes * spatial_dimension,
                  int_t_N_val + size_of_shapes * spatial_dimension);

      for (Int n = 0; n < size_of_shapes * spatial_dimension; ++n) {
        int_t_N_val[n] *= -1.;
      }

      int_t_N_val += nb_nodes_per_element * spatial_dimension;
      partial_int_t_N_val += size_of_shapes * spatial_dimension;
    }

    /**
     * compute @f$\int err \cdot N\, dS@f$ by  @f$ \sum_q \mathbf{N}^t
     * \mathbf{err}_q \overline w_q J_q@f$
     */
    auto int_err_N = std::make_shared<Array<Real>>(
        nb_element, spatial_dimension * size_of_shapes, "int_err_N");

    fem_cohesive.integrate(*err_opening_cpy, *int_err_N,
                           spatial_dimension * size_of_shapes, type, ghost_type,
                           elem_filter);

    /// assemble
    model->getDOFManager().assembleElementalArrayLocalArray(
        "displacement", *int_t_N, internal_force, type, ghost_type, 1,
        elem_filter);

    //    auto lambda_connectivity = lambda_connectivities(type, ghost_type);
    auto underlying_type = Mesh::getFacetType(type);
    model->getDOFManager().assembleElementalArrayToResidual(
        "lambda", *int_err_N, underlying_type, ghost_type, 1., elem_filter);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::assembleStiffnessMatrix(
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // auto & lambda_connectivities =
  //     model->getMesh().getElementalData<Idx>("lambda_connectivities");

  for (auto type : getElementFilter().elementTypes(spatial_dimension,
                                                   ghost_type, _ek_cohesive)) {
    auto nb_quadrature_points =
        fem_cohesive.getNbIntegrationPoints(type, ghost_type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    const auto & elem_filter = getElementFilter(type, ghost_type);
    auto nb_element = elem_filter.size();

    if (nb_element == 0U) {
      continue;
    }

    const auto & shapes = fem_cohesive.getShapes(type, ghost_type);
    auto size_of_shapes = shapes.getNbComponent();

    auto shapes_filtered = std::make_shared<Array<Real>>(
        nb_element * nb_quadrature_points, size_of_shapes, "filtered shapes");

    for (auto && data :
         zip(filter(elem_filter,
                    make_view(shapes, size_of_shapes, nb_quadrature_points)),
             make_view(*shapes_filtered, size_of_shapes,
                       nb_quadrature_points))) {
      std::get<1>(data) = std::get<0>(data);
    }

    Matrix<Real> A(spatial_dimension * size_of_shapes,
                   spatial_dimension * nb_nodes_per_element);
    A.zero();
    for (Int i = 0; i < spatial_dimension * size_of_shapes; ++i) {
      A(i, i) = 1;
      A(i, i + spatial_dimension * size_of_shapes) = -1;
    }

    /// get the tangent matrix @f$\frac{\partial{(t/\delta)}}{\partial{\delta}}
    /// @f$
    /// TODO : optimisation not to reassemble uu term, which does not change
    /// during computation
    auto tangent_stiffness_matrix_uu = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        spatial_dimension * spatial_dimension, "tangent_stiffness_matrix_uu");

    auto tangent_stiffness_matrix_ll = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        spatial_dimension * spatial_dimension, "tangent_stiffness_matrix_ll");

    computeNormal(model->getCurrentPosition(), normals(type, ghost_type), type,
                  ghost_type);

    /// compute openings @f$\mathbf{\delta}@f$
    // computeOpening(model->getDisplacement(), opening(type, ghost_type), type,
    // ghost_type);

    tangent_stiffness_matrix_uu->zero();
    tangent_stiffness_matrix_ll->zero();

    computeTangentTraction(type, *tangent_stiffness_matrix_uu,
                           *tangent_stiffness_matrix_ll, ghost_type);

    UInt size_at_nt_duu_n_a = spatial_dimension * nb_nodes_per_element *
                              spatial_dimension * nb_nodes_per_element;
    auto at_nt_duu_n_a =
        std::make_unique<Array<Real>>(nb_element * nb_quadrature_points,
                                      size_at_nt_duu_n_a, "A^t*N^t*Duu*N*A");

    UInt size_nt_dll_n =
        spatial_dimension * size_of_shapes * spatial_dimension * size_of_shapes;
    auto nt_dll_n = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, size_nt_dll_n, "N^t*Dll*N");

    UInt size_at_nt_dul_n = spatial_dimension * nb_nodes_per_element *
                            spatial_dimension * size_of_shapes;
    auto at_nt_dul_n = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, size_at_nt_dul_n, "A^t*N^t*Dul*N");

    Matrix<Real> N(spatial_dimension, spatial_dimension * size_of_shapes);

    for (auto && [At_Nt_Duu_N_A, Duu, Nt_Dll_N, Dll, At_Nt_Dul_N, shapes] :
         zip(make_view(*at_nt_duu_n_a, spatial_dimension * nb_nodes_per_element,
                       spatial_dimension * nb_nodes_per_element),
             make_view(*tangent_stiffness_matrix_uu, spatial_dimension,
                       spatial_dimension),
             make_view(*nt_dll_n, spatial_dimension * size_of_shapes,
                       spatial_dimension * size_of_shapes),
             make_view(*tangent_stiffness_matrix_ll, spatial_dimension,
                       spatial_dimension),
             make_view(*at_nt_dul_n, spatial_dimension * nb_nodes_per_element,
                       spatial_dimension * size_of_shapes),
             make_view(*shapes_filtered, size_of_shapes))) {
      N.zero();
      /**
       * store  the   shapes  in  voigt   notations  matrix  @f$\mathbf{N}  =
       * \begin{array}{cccccc} N_0(\xi) & 0 & N_1(\xi)  &0 & N_2(\xi) & 0 \\
       * 0 & * N_0(\xi)& 0 &N_1(\xi)& 0 & N_2(\xi) \end{array} @f$
       **/
      for (Int i = 0; i < spatial_dimension; ++i) {
        for (Int n = 0; n < size_of_shapes; ++n) {
          N(i, i + spatial_dimension * n) = shapes(n);
        }
      }

      /**
       * compute stiffness matrix  @f$   \mathbf{K}    =    \delta
       *\mathbf{U}^T \int_{\Gamma_c}    {\mathbf{P}^t
       *\frac{\partial{\mathbf{t}}}
       *{\partial{\delta}}
       * \mathbf{P} d\Gamma \Delta \mathbf{U}}  @f$
       **/
      auto && NA = N * A;
      At_Nt_Duu_N_A = (Duu * NA).transpose() * NA;
      Nt_Dll_N = (Dll * N).transpose() * N;
      At_Nt_Dul_N = NA.transpose() * N;
    }

    auto Kuu_e =
        std::make_unique<Array<Real>>(nb_element, size_at_nt_duu_n_a, "Kuu_e");

    fem_cohesive.integrate(*at_nt_duu_n_a, *Kuu_e, size_at_nt_duu_n_a, type,
                           ghost_type, elem_filter);

    auto Kll_e =
        std::make_unique<Array<Real>>(nb_element, size_nt_dll_n, "Kll_e");

    fem_cohesive.integrate(*nt_dll_n, *Kll_e, size_nt_dll_n, type, ghost_type,
                           elem_filter);

    auto Kul_e =
        std::make_unique<Array<Real>>(nb_element, size_at_nt_dul_n, "Kul_e");

    fem_cohesive.integrate(*at_nt_dul_n, *Kul_e, size_at_nt_dul_n, type,
                           ghost_type, elem_filter);

    model->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "displacement", *Kuu_e, type, ghost_type, _symmetric, elem_filter);

    auto underlying_type = Mesh::getFacetType(type);
    model->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "lambda", *Kll_e, underlying_type, ghost_type, _symmetric,
        elem_filter);

    auto connectivity = model->getMesh().getConnectivity(type, ghost_type);
    auto conn = make_view(connectivity, connectivity.getNbComponent()).begin();

    auto lambda_connectivity =
        model->getLambdaMesh().getConnectivity(underlying_type, ghost_type);
    auto lambda_conn =
        make_view(lambda_connectivity, lambda_connectivity.getNbComponent())
            .begin();

    /// Assemble Kll_e
    // TermsToAssemble term_ll("lambda", "lambda");
    // auto el_mat_it_ll = Kll_e->begin(spatial_dimension * size_of_shapes,
    //                                  spatial_dimension * size_of_shapes);

    // auto compute_ll = [&](const auto & el) {
    //   auto kll_e = *el_mat_it_ll;
    //   auto && lda_conn_el = lambda_conn[el];
    //   auto N = lda_conn_el.rows();
    //   for (Int m = 0; m < N; ++m) {
    //     auto ldai = lda_conn_el(m);
    //     for (Int n = m; n < N; ++n) {
    //       auto ldaj = lda_conn_el(n);
    //       for (Int k = 0; k < spatial_dimension; ++k) {
    //         for (Int l = 0; l < spatial_dimension; ++l) {
    //           auto && term_ll_ij = term_ll(ldai * spatial_dimension + k,
    //                                        ldaj * spatial_dimension + l);
    //           term_ll_ij =
    //               kll_e(m * spatial_dimension + k, n * spatial_dimension +
    //               l);
    //         }
    //       }
    //     }
    //   }
    //   ++el_mat_it_ll;
    // };
    // for_each_element(nb_element, elem_filter, compute_ll);

    // model->getDOFManager().assemblePreassembledMatrix("K", term_ll);
    // model->getDOFManager().getMatrix("K").saveMatrix("Kuu_terms.mtx");

    /// Assemble Klu_e
    TermsToAssemble term_ul("displacement", "lambda");
    auto el_mat_it_ul = Kul_e->begin(spatial_dimension * nb_nodes_per_element,
                                     spatial_dimension * size_of_shapes);

    auto compute_ul = [&](const auto & el) {
      auto kul_e = *el_mat_it_ul;
      auto && u_conn_el = conn[el];
      auto && lda_conn_el = lambda_conn[el];
      auto M = u_conn_el.rows();
      auto N = lda_conn_el.rows();
      for (Int m = 0; m < M; ++m) {
        for (Int n = 0; n < N; ++n) {
          auto u = u_conn_el(m);
          auto lda = lda_conn_el(n);
          for (Int k = 0; k < spatial_dimension; ++k) {
            for (Int l = 0; l < spatial_dimension; ++l) {
              auto && term_ul_ij = term_ul(u * spatial_dimension + k,
                                           lda * spatial_dimension + l);
              term_ul_ij =
                  kul_e(m * spatial_dimension + k, n * spatial_dimension + l);
            }
          }
        }
      }
      ++el_mat_it_ul;
    };
    for_each_element(nb_element, elem_filter, compute_ul);

    model->getDOFManager().assemblePreassembledMatrix("K", term_ul);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::computeLambdaOnQuad(ElementType type,
                                                      GhostType ghost_type) {
  auto & fem_lambda = this->model->getFEEngine("LambdaFEEngine");
  const auto & lambda = this->model->getLambda();
  auto & lambda_on_quad = this->lambda(type, ghost_type);

  std::cout << "lambda in : MaterialCohesiveDamage<dim>::computeLambdaOnQuad " << std::endl;
  lambda.printself(std::cout << std::endl);
  ArrayPrintHelper<true>::print_content(lambda,std::cout,0);


  auto underlying_type = Mesh::getFacetType(type);
  fem_lambda.interpolateOnIntegrationPoints(
      lambda, lambda_on_quad, dim, underlying_type, ghost_type,
      this->getElementFilter(type, ghost_type));

  std::cout << "lambda_on_quad in : MaterialCohesiveDamage<dim>::computeLambdaOnQuad " << std::endl;
  lambda_on_quad.printself(std::cout << std::endl);
  ArrayPrintHelper<true>::print_content(lambda_on_quad,std::cout,0);

  std::cout << "shapes_lambda in : MaterialCohesiveDamage<dim>::computeLambdaOnQuad " << std::endl;
  auto shapes_lambda = fem_lambda.getShapes(underlying_type,ghost_type,0);
  ArrayPrintHelper<true>::print_content(shapes_lambda,std::cout,0);

  std::cout << "shapes_cohesive in : MaterialCohesiveDamage<dim>::computeLambdaOnQuad " << std::endl;
  auto & fem_cohesive =
      this->model->getFEEngineClass<MyFEEngineCohesiveType>("CohesiveFEEngine");
  auto shapes_cohesive = fem_cohesive.getShapes(type,ghost_type,0);
  ArrayPrintHelper<true>::print_content(shapes_cohesive,std::cout,0);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::computeTraction(ElementType el_type,
                                                  GhostType ghost_type) {

  for (const auto & type : getElementFilter().elementTypes(
           spatial_dimension, ghost_type, _ek_cohesive)) {
    auto & elem_filter = getElementFilter(type, ghost_type);
    auto nb_element = elem_filter.size();
    if (nb_element == 0) {
      continue;
    }

    /// compute normals @f$\mathbf{n}@f$
    computeNormal(model->getCurrentPosition(), normals(type, ghost_type), type,
                  ghost_type);

    /// compute openings @f$\mathbf{\delta}@f$
    computeOpening(model->getDisplacement(), opening(type, ghost_type), type,
                   ghost_type);

    computeLambdaOnQuad(el_type, ghost_type);

    auto & traction = tractions(type, ghost_type);
    std::cout << "traction in : MaterialCohesiveDamage<dim>::computeTraction " << std::endl;
    ArrayPrintHelper<true>::print_content(traction,std::cout,0);

    auto & lambda_ = lambda(type, ghost_type);
    std::cout << "lambda in : MaterialCohesiveDamage<dim>::computeTraction " << std::endl;
    ArrayPrintHelper<true>::print_content(lambda_,std::cout,0);

    for (auto && args : getArguments(el_type, ghost_type)) {
      this->computeTractionOnQuad(args);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::computeTangentTraction(
    ElementType el_type, Array<Real> & tangent_matrix_uu,
    Array<Real> & tangent_matrix_ll, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && [args, tangent_uu, tangent_ll] :
       zip(getArguments(el_type, ghost_type),
           make_view<dim, dim>(tangent_matrix_uu),
           make_view<dim, dim>(tangent_matrix_ll))) {
    computeTangentTractionOnQuad(tangent_uu, tangent_ll, args);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template class MaterialCohesiveDamage<1>;
template class MaterialCohesiveDamage<2>;
template class MaterialCohesiveDamage<3>;
const bool material_is_alocated_cohesive_damage [[maybe_unused]] =
    instantiateMaterial<MaterialCohesiveDamage>("cohesive_damage");

} // namespace akantu
