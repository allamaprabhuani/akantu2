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
    : MaterialCohesive(model, id), lambda("lambda",*this) {
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
  lambda.initialize(dim);

  const auto & mesh_facets = model->getMeshFacets();
  for (const auto & type_facet : mesh_facets.elementTypes(dim - 1)) {
    auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);

    const auto & facet_filter_array = facet_filter(type_facet);
    const auto & lambda_array = lambda(type_cohesive);

    auto nb_quad_facet =
        model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);

    for (auto && [facet, lda] :
         zip(facet_filter_array, lambda_array)) {

//        for (Int q = 0; q < nb_quad_facet; ++q) {
//            auto current_quad = facet * nb_quad_facet + q;
//        }
        std::cout << " lda = " << lda << std::endl;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamage<dim>::computeTraction(ElementType el_type,
                                                  GhostType ghost_type) {
  std::cout << "Print MaterialCohesiveDamage<dim>::computeTraction : TODO " << std::endl;
  throw;
}

/* -------------------------------------------------------------------------- */
template class MaterialCohesiveDamage<1>;
template class MaterialCohesiveDamage<2>;
template class MaterialCohesiveDamage<3>;
static bool material_is_alocated_cohesive_damage =
    instantiateMaterial<MaterialCohesiveDamage>("cohesive_damage");

} // namespace akantu
