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
      lambda(registerInternal<Real, CohesiveInternalField>("lambda", spatial_dimension)),
      err_openings(registerInternal<Real, CohesiveInternalField>("err_openings", spatial_dimension)),
      czm_damage(registerInternal<Real, CohesiveInternalField>("czm_damage", 1))
      {
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
  auto & lambda_connectivities =
          model->getMesh().getElementalData<Idx>("lambda_connectivities");

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

      for (Int q = 0; q < nb_quadrature_points; ++q,
           ++traction_it,++err_opening_it,
           ++current_quad,
           ++traction_cpy_it,++err_opening_cpy_it) {

        auto && shapes_filtered = shapes_filtered_begin[current_quad];

        *traction_cpy_it =
            (*traction_it ) * shapes_filtered;
        *err_opening_cpy_it =
            (*err_opening_it ) * shapes_filtered;
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

    /// TODO : not sure this should be done using fem_cohesive ?
    fem_cohesive.integrate(*err_opening_cpy, *int_err_N,
                           spatial_dimension * size_of_shapes, type, ghost_type,
                           elem_filter);

    /// assemble
    model->getDOFManager().assembleElementalArrayLocalArray(
        *int_t_N, internal_force, type, ghost_type, 1, elem_filter);

    /// Not sure why we need connectivity here ?
    auto lambda_connectivity = lambda_connectivities(type, ghost_type);
    model->getDOFManager().assembleElementalArrayToResidual("lambda",*int_err_N,
                                                           lambda_connectivity,
                                                           type,ghost_type,1.,elem_filter);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & lambda_connectivities =
      model->getMesh().getElementalData<Idx>("lambda_connectivities");

  for (auto type : getElementFilter().elementTypes(spatial_dimension,
                                                   ghost_type, _ek_cohesive)) {
    auto nb_quadrature_points =
        fem_cohesive.getNbIntegrationPoints(type, ghost_type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    auto & elem_filter = getElementFilter(type, ghost_type);
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

    computeTangentTraction(type, *tangent_stiffness_matrix_uu, *tangent_stiffness_matrix_ll, ghost_type);

    UInt size_at_nt_duu_n_a = spatial_dimension * nb_nodes_per_element *
                            spatial_dimension * nb_nodes_per_element;
    auto at_nt_duu_n_a = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, size_at_nt_duu_n_a, "A^t*N^t*Duu*N*A");

    UInt size_nt_dll_n = spatial_dimension * nb_nodes_per_element *
                         spatial_dimension * nb_nodes_per_element;
    auto nt_dll_n = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, size_nt_dll_n, "N^t*Dll*N");

    UInt size_at_nt_dul_n = spatial_dimension * nb_nodes_per_element *
                            spatial_dimension * nb_nodes_per_element;
    auto at_nt_dul_n = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, size_at_nt_dul_n, "A^t*N^t*Dul*N");

    Matrix<Real> N(spatial_dimension, spatial_dimension * nb_nodes_per_element);

    for (auto && data :
         zip(make_view(*at_nt_duu_n_a, spatial_dimension * nb_nodes_per_element,
                       spatial_dimension * nb_nodes_per_element),
             make_view(*tangent_stiffness_matrix_uu, spatial_dimension,
                       spatial_dimension),
             make_view(*nt_dll_n, spatial_dimension * nb_nodes_per_element,
                       spatial_dimension * nb_nodes_per_element),
             make_view(*tangent_stiffness_matrix_ll, spatial_dimension,
                       spatial_dimension),
             make_view(*at_nt_dul_n, spatial_dimension * nb_nodes_per_element,
                       spatial_dimension * nb_nodes_per_element),
             make_view(*shapes_filtered, size_of_shapes))) {

      auto && At_Nt_Duu_N_A = std::get<0>(data);
      auto && Duu = std::get<1>(data);
      auto && Nt_Dll_N = std::get<2>(data);
      auto && Dll = std::get<3>(data);
      auto && At_Nt_Dul_N = std::get<4>(data);
      auto && shapes = std::get<5>(data);
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
      At_Nt_Dul_N =  NA.transpose() * N;
    }

    auto Kuu_e =
        std::make_unique<Array<Real>>(nb_element, size_at_nt_duu_n_a, "Kuu_e");

    fem_cohesive.integrate(*at_nt_duu_n_a, *Kuu_e, size_at_nt_duu_n_a, type,
                           ghost_type, elem_filter);

    auto Kll_e =
        std::make_unique<Array<Real>>(nb_element, size_nt_dll_n, "Kll_e");

    fem_cohesive.integrate(*nt_dll_n, *Kll_e, size_nt_dll_n, type,
                           ghost_type, elem_filter);

    auto Kul_e =
        std::make_unique<Array<Real>>(nb_element, size_at_nt_dul_n, "Kul_e");

    fem_cohesive.integrate(*at_nt_dul_n, *Kul_e, size_at_nt_dul_n, type,
                           ghost_type, elem_filter);

    model->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "displacement", *Kuu_e, type, ghost_type, _symmetric, elem_filter);


    /// Not sure why we need connectivity here ?
    auto lambda_connectivity = lambda_connectivities(type, ghost_type);
    model->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "lambda", *Kll_e, lambda_connectivity,type, ghost_type, _symmetric, elem_filter);


    /// Where do we tell TermsToAssemble to use Kul_e ???
    /// Do we need to assemble Klu_e
    TermsToAssemble term_ul("displacement","lambda");

    auto conn_it = lambda_connectivity.begin(nb_nodes_per_element);
    auto el_mat_it = Kul_e->begin(spatial_dimension * nb_nodes_per_element,
                                  spatial_dimension * nb_nodes_per_element);


    /// TODO : using lambda connectivity to properly assemble ul terms
    for (Int el = 0; el < nb_element; ++el, ++conn_it, ++el_mat_it) {
        std::cout << "conn_it = " << *conn_it << std::endl;
    }
//    for(Int i = ..., Int j = ...)
//    {
//        TermToAssemble term_ul_ij(i,j);
//        term_ul_ij = ...
//        term_ul(i,j) =  term_ul_ij;
//    }
    model->getDOFManager().assemblePreassembledMatrix("K",term_ul);
  }

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
  ///         WILL REQUIRE LAMBDA CONNECTIVITY
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
template <Int dim>
void MaterialCohesiveDamage<dim>::computeTangentTraction(ElementType el_type,
                                                         Array<Real> & tangent_matrix_uu,
                                                         Array<Real> &tangent_matrix_ll,
                                                         GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && [args, tangent_uu, tangent_ll] : zip(getArguments(el_type, ghost_type),
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
