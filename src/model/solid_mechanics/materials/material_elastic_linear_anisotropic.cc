/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "material_elastic_linear_anisotropic.hh"
#include "solid_mechanics_model.hh"
#include <algorithm>
#include <sstream>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialElasticLinearAnisotropic<dim>::MaterialElasticLinearAnisotropic(
    SolidMechanicsModel & model, const ID & id, bool symmetric)
    : Material(model, id), rot_mat(dim, dim), Cprime(dim * dim, dim * dim),
      C(voigt_h::size, voigt_h::size), eigC(voigt_h::size),
      symmetric(symmetric), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();

  for (int i : arange(dim)) {
    this->dir_vecs.emplace_back(std::make_unique<Vector<Real, dim>>());
    auto && n = *this->dir_vecs.back();
    n.zero();
    n[i] = 1.;
    this->registerParam("n" + std::to_string(i + 1), *(this->dir_vecs.back()),
                        _pat_parsmod, "Direction of main material axis");
  }

  for (auto i : arange(voigt_h::size)) {
    decltype(i) start = 0;
    if (this->symmetric) {
      start = i;
    }
    for (auto j : arange(start, voigt_h::size)) {
      auto param = "C" + std::to_string(i + 1) + std::to_string(j + 1);
      this->registerParam(param, this->Cprime(i, j), Real(0.), _pat_parsmod,
                          "Coefficient " + param);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialElasticLinearAnisotropic<dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElasticLinearAnisotropic<dim>::updateInternalParameters() {

  Material::updateInternalParameters();
  if (this->symmetric) {
    for (auto i : arange(voigt_h::size)) {
      for (auto j : arange(i + 1, voigt_h::size)) {
        this->Cprime(j, i) = this->Cprime(i, j);
      }
    }
  }
  this->rotateCprime();
  this->C.eig(this->eigC);

  this->was_stiffness_assembled = false;
}

/* -------------------------------------------------------------------------- */
template <Int Dim> void MaterialElasticLinearAnisotropic<Dim>::rotateCprime() {
  // start by filling the empty parts fo Cprime
  auto diff = Dim * Dim - voigt_h::size;
  for (auto i : arange(voigt_h::size, Dim * Dim)) {
    for (auto j : arange(Dim * Dim)) {
      this->Cprime(i, j) = this->Cprime(i - diff, j);
    }
  }
  for (auto i : arange(Dim * Dim)) {
    for (auto j : arange(voigt_h::size, Dim * Dim)) {
      this->Cprime(i, j) = this->Cprime(i, j - diff);
    }
  }
  // construction of rotator tensor
  // normalise rotation matrix
  for (auto j : arange(Dim)) {
    auto && rot_vec = this->rot_mat(j);
    rot_vec = *this->dir_vecs[j];
    rot_vec.normalize();
  }

  // make sure the vectors form a right-handed base
  Real test_axis = 0.;
  if (Dim == 2) {
    Vector<Real, 3> v1 = Vector<Real, 3>::Zero();
    Vector<Real, 3> v2 = Vector<Real, 3>::Zero();
    v1.block<Dim, 1>(0, 0) = this->rot_mat(0);
    v2.block<Dim, 1>(0, 0) = this->rot_mat(1);

    Vector<Real, 3> v3 = v1.cross(v2);
    if (v3.norm() < 8 * std::numeric_limits<Real>::epsilon()) {
      AKANTU_ERROR("The axis vectors parallel.");
    }

    v3.normalize();
    test_axis = (v1.cross(v2) - v3).norm();
  } else if (Dim == 3) {
    Vector<Real, 3> v1 = this->rot_mat(0);
    Vector<Real, 3> v2 = this->rot_mat(1);
    Vector<Real, 3> v3 = this->rot_mat(2);
    test_axis = (v1.cross(v2) - v3).norm();
  }

  if (test_axis > 8 * std::numeric_limits<Real>::epsilon()) {
    AKANTU_ERROR("The axis vectors do not form a right-handed coordinate "
                 << "system. I. e., ||n1 x n2 - n3|| should be zero, but "
                 << "it is " << test_axis << ".");
  }

  // create the rotator and the reverse rotator
  Matrix<Real, Dim * Dim, Dim * Dim> rotator;
  Matrix<Real, Dim * Dim, Dim * Dim> revrotor;
  for (auto i : arange(Dim)) {
    for (auto j : arange(Dim)) {
      for (auto k : arange(Dim)) {
        for (auto l : arange(Dim)) {
          auto I = voigt_h::mat[i][j];
          auto J = voigt_h::mat[k][l];
          rotator(I, J) = this->rot_mat(k, i) * this->rot_mat(l, j);
          revrotor(I, J) = this->rot_mat(i, k) * this->rot_mat(j, l);
        }
      }
    }
  }

  // create the full rotated matrix
  Matrix<Real, Dim * Dim, Dim * Dim> Cfull(Dim * Dim, Dim * Dim);
  Cfull = rotator * Cprime * revrotor;

  for (auto i : arange(voigt_h::size)) {
    for (auto j : arange(voigt_h::size)) {
      this->C(i, j) = Cfull(i, j);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElasticLinearAnisotropic<dim>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  auto && arguments = getArguments(el_type, ghost_type);
  for (auto && data : arguments) {
    this->computeStressOnQuad(data);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElasticLinearAnisotropic<dim>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  auto && arguments =
      Material::getArgumentsTangent<dim>(tangent_matrix, el_type, ghost_type);

  for (auto && args : arguments) {
    this->computeTangentModuliOnQuad(args);
  }

  this->was_stiffness_assembled = true;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElasticLinearAnisotropic<dim>::computePotentialEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_ASSERT(!this->finite_deformation,
                      "finite deformation not possible in material anisotropic "
                      "(TO BE IMPLEMENTED)");

  auto && arguments = Material::getArguments<dim>(el_type, _not_ghost);

  for (auto && args :
       zip(arguments, this->potential_energy(el_type, _not_ghost))) {
    this->computePotentialEnergyOnQuad(std::get<0>(args), std::get<1>(args));
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialElasticLinearAnisotropic<dim>::computePotentialEnergyByElement(
    ElementType type, Int index, Vector<Real> & epot_on_quad_points) {

  auto gradu_view = make_view<dim, dim>(this->gradu(type));
  auto stress_view = make_view<dim, dim>(this->stress(type));

  auto nb_quadrature_points = this->fem.getNbIntegrationPoints(type);

  auto gradu_it = gradu_view.begin() + index * nb_quadrature_points;
  auto gradu_end = gradu_it + nb_quadrature_points;
  auto stress_it = stress_view.begin() + index * nb_quadrature_points;
  auto stress_end = stress_it + nb_quadrature_points;

  auto epot_quad = epot_on_quad_points.begin();

  Matrix<Real> grad_u(dim, dim);

  for (auto data : zip("grad_u"_n = range(gradu_it, gradu_end),
                   "sigma"_n = range(stress_it, stress_end),
                   "Epot"_n = epot_on_quad_points)) {
    this->computePotentialEnergyOnQuad(data, data["Epot"_n]);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialElasticLinearAnisotropic<dim>::getCelerity(
    const Element & /*element*/) const {
  return std::sqrt(this->eigC(0) / rho);
}

/* -------------------------------------------------------------------------- */
template class MaterialElasticLinearAnisotropic<1>;
template class MaterialElasticLinearAnisotropic<2>;
template class MaterialElasticLinearAnisotropic<3>;

static bool material_is_alocated_elastic =
    instantiateMaterial<MaterialElasticLinearAnisotropic>(
        "elastic_anisotropic");

} // namespace akantu
