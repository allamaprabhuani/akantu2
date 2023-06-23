/**
 * @file   material_bulk_viscosity.hh
 *
 * @author Shenghan Zhang <shenghan.zhang@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Jul 20 2016
 *
 * @brief  Bulk viscosity in the material
 *
 * @section LICENSE
 *
 * Copyright (©) 2016 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory EESD (Earthquake Engineering and Structural Dynamics)
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
#include "material_bulk_viscosity.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, class Parent>
MaterialBulkViscosity<dim, Parent>::MaterialBulkViscosity(
    SolidMechanicsModel & model, const ID & id)
    : Parent(model, id), characteristic_length("characteristic_length", *this),
      sound_speed("sound_speed", *this), grad_v("grad_v", *this),
      sigma_elastic("sigma_elastic", *this),
      sigma_viscous("sigma_viscous", *this) {
  this->registerParam("b1", b1, 0.06, _pat_parsable | _pat_readable,
                      "Linear damping coefficient");
  this->registerParam("b2", b2, 1.2, _pat_parsable | _pat_readable,
                      "Quadratic damping coefficient");

  this->registerParam("apply_always", apply_always, true,
                      _pat_parsable | _pat_readable,
                      "Applied always or only on negative strain rates");
}

/* -------------------------------------------------------------------------- */
template <Int dim, class Parent>
void MaterialBulkViscosity<dim, Parent>::initMaterial() {
  Parent::initMaterial();

  this->characteristic_length.initialize(1);
  this->sound_speed.initialize(1);
  this->grad_v.initialize(dim * dim);

  this->sigma_elastic.initialize(dim * dim);
  this->sigma_viscous.initialize(dim * dim);

  this->updateInternals();
}

/* -------------------------------------------------------------------------- */
template <Int dim, class Parent>
void MaterialBulkViscosity<dim, Parent>::updateInternals() {
  Mesh & mesh = this->fem.getMesh();

  for (auto ghost_type : ghost_types) {
    for (auto type : this->element_filter.elementTypes(dim, ghost_type)) {
      const auto & elem_filter = this->element_filter(type, ghost_type);
      auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);

      Element elem{type, 0, ghost_type};

      Array<Real> Xs(0, nb_nodes_per_element * dim);
      FEEngine::extractNodalToElementField(mesh,
                                           this->model.getCurrentPosition(), Xs,
                                           type, ghost_type, elem_filter);

      for (auto && [X, Le, c] :
           zip(make_view(Xs, dim, nb_nodes_per_element),
               this->characteristic_length(type, ghost_type),
               this->sound_speed(type, ghost_type))) {
        Le = this->fem.getElementInradius(X, type);
        c = this->getCelerity(elem);
        ++elem.element;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, class Parent>
template <typename Args>
inline void
MaterialBulkViscosity<dim, Parent>::computeStressOnQuad(Args && args) const {
  auto && grad_v = args["grad_v"_n];
  auto && sigma = args["sigma"_n];
  auto && le = args["l_e"_n];
  auto && c = args["sound_speed"_n];

  Real tr_D = grad_v.trace();

  if (tr_D < 0. or apply_always) {
    // Quadratic contribution to spread out the wave front to account
    // for shock discontinuities
    Real p_bv2 = 0.;
    if (tr_D < 0.) {
      p_bv2 = this->rho * Math::pow<2>(this->b2 * le * tr_D);
    }
    // Linear contribution to account for the post-shocks instabilities
    Real p_bv1 = this->rho * this->b1 * le * c * tr_D;

    Real p_bv = p_bv2 - p_bv1;
    sigma.eye(p_bv);
  } else {
    sigma.clear();
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, class Parent>
void MaterialBulkViscosity<dim, Parent>::computeStress(ElementType el_type,
                                                       GhostType ghost_type) {
  Parent::computeStress(el_type, ghost_type);

  this->fem.gradientOnIntegrationPoints(
      this->model.getVelocity(), this->grad_v(el_type, ghost_type), dim,
      el_type, ghost_type, this->element_filter(el_type, ghost_type));

  auto && arguments = getArguments(el_type, ghost_type);
  for (auto && args : arguments) {
    this->computeStressOnQuad(args);
    args["sigma_elastic"_n] = args["sigma"_n];
    args["sigma"_n] -= args["sigma_viscous"_n];
  }
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
