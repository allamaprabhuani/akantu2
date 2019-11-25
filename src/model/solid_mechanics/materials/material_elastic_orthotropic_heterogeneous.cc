/**
 * @file   material_elastic_orthotropic.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Orthotropic elastic material
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "material_elastic_orthotropic_heterogeneous.hh"
#include "solid_mechanics_model.hh"
#include <algorithm>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt Dim>
MaterialElasticOrthotropicHeterogeneous<
    Dim>::MaterialElasticOrthotropicHeterogeneous(SolidMechanicsModel & model,
                                                  const ID & id)
    : MaterialElasticLinearAnisotropicHeterogeneous<Dim>(model, id),
      E1_field("E1_field", *this), E2_field("E2_field", *this),
      E3_field("E3_field", *this), nu12_field("nu12_field", *this),
      nu13_field("nu13_field", *this), nu23_field("nu23_field", *this),
      G12_field("G12_field", *this), G13_field("G13_field", *this),
      G23_field("G23_field", *this), are_internals_initialized(false) {
  AKANTU_DEBUG_IN();
  this->registerParam("E1", E1, Real(0.), _pat_parsmod, "Young's modulus (n1)");
  this->registerParam("E2", E2, Real(0.), _pat_parsmod, "Young's modulus (n2)");
  this->registerParam("nu12", nu12, Real(0.), _pat_parsmod,
                      "Poisson's ratio (12)");
  this->registerParam("G12", G12, Real(0.), _pat_parsmod, "Shear modulus (12)");

  this->registerParam("E3", E3, Real(0.), _pat_parsmod, "Young's modulus (n3)");
  this->registerParam("nu13", nu13, Real(0.), _pat_parsmod,
                      "Poisson's ratio (13)");
  this->registerParam("nu23", nu23, Real(0.), _pat_parsmod,
                      "Poisson's ratio (23)");
  this->registerParam("G13", G13, Real(0.), _pat_parsmod, "Shear modulus (13)");
  this->registerParam("G23", G23, Real(0.), _pat_parsmod, "Shear modulus (23)");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt Dim>
void MaterialElasticOrthotropicHeterogeneous<Dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  parent::initMaterial();

  this->E1_field.setDefaultValue(this->E1);
  this->E2_field.setDefaultValue(this->E2);
  this->E3_field.setDefaultValue(this->E3);
  this->nu12_field.setDefaultValue(this->nu12);
  this->nu13_field.setDefaultValue(this->nu13);
  this->nu23_field.setDefaultValue(this->nu23);
  this->G12_field.setDefaultValue(this->G12);
  this->G13_field.setDefaultValue(this->G13);
  this->G23_field.setDefaultValue(this->G23);
  this->E1_field.initialize(1);
  this->E2_field.initialize(1);
  this->nu12_field.initialize(1);
  this->G12_field.initialize(1);
  this->E3_field.initialize(1);
  this->nu13_field.initialize(1);
  this->nu23_field.initialize(1);
  this->G13_field.initialize(1);
  this->G23_field.initialize(1);

  this->are_internals_initialized = true;
  AKANTU_DEBUG_ASSERT(not this->finite_deformation,
                      "finite deformation not possible in material orthotropic "
                      "(TO BE IMPLEMENTED)");

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt Dim>
void MaterialElasticOrthotropicHeterogeneous<Dim>::updateInternalParameters() {
  if (not this->are_internals_initialized)
    return;

  /* 1) construction of temporary material frame stiffness tensor------------ */
  // http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {

    auto & E1_f = this->E1_field(el_type, _not_ghost);
    auto & E2_f = this->E2_field(el_type, _not_ghost);
    auto & E3_f = this->E3_field(el_type, _not_ghost);
    auto & nu12_f = this->nu12_field(el_type, _not_ghost);
    auto & nu13_f = this->nu13_field(el_type, _not_ghost);
    auto & nu23_f = this->nu23_field(el_type, _not_ghost);
    auto & G12_f = this->G12_field(el_type, _not_ghost);
    auto & G13_f = this->G13_field(el_type, _not_ghost);
    auto & G23_f = this->G23_field(el_type, _not_ghost);
    auto & Cprime_f = this->Cprime_field(el_type, _not_ghost);
    auto & C_f = this->C_field(el_type, _not_ghost);
    auto & eigC_f = this->eigC_field(el_type, _not_ghost);
    auto & dir_vecs_f = this->dir_vecs_field(el_type, _not_ghost);

    for (auto && data :
         zip(E1_f, E2_f, E3_f, nu12_f, nu13_f, nu23_f, G12_f, G13_f, G23_f,
             make_view(Cprime_f, Dim * Dim, Dim * Dim),
             make_view(C_f, voigt_h::size, voigt_h::size),
             make_view(eigC_f, voigt_h::size),
             make_view(dir_vecs_f, Dim, Dim))) {

      auto & _E1 = std::get<0>(data);
      auto & _E2 = std::get<1>(data);
      auto & _E3 = std::get<2>(data);
      auto & _nu12 = std::get<3>(data);
      auto & _nu13 = std::get<4>(data);
      auto & _nu23 = std::get<5>(data);
      auto & _G12 = std::get<6>(data);
      auto & _G13 = std::get<7>(data);
      auto & _G23 = std::get<8>(data);
      auto & _Cprime = std::get<9>(data);
      auto & _C = std::get<10>(data);
      auto & _eigC = std::get<11>(data);
      auto & _dir_vecs = std::get<12>(data);

      updateInternalParametersOnQuad(_E1, _E2, _E3, _nu12, _nu13, _nu23, _G12,
                                     _G13, _G23, _Cprime, _C, _eigC, _dir_vecs);
    }
  }
}
/* --------------------------------------------------------------------------
 */
template <UInt Dim>
void MaterialElasticOrthotropicHeterogeneous<
    Dim>::updateInternalParametersOnQuad(const Real & _E1, const Real & _E2,
                                         const Real & _E3, const Real & _nu12,
                                         const Real & _nu13, const Real & _nu23,
                                         const Real & _G12, const Real & _G13,
                                         const Real & _G23,
                                         Matrix<Real> & _Cprime,
                                         Matrix<Real> & _C,
                                         Vector<Real> & _eigC,
                                         Matrix<Real> & _dir_vecs) {

  Real _nu21 = _nu12 * _E2 / _E1;
  Real _nu31 = _nu13 * _E3 / _E1;
  Real _nu32 = _nu23 * _E3 / _E2;

  // Full (i.e. dim^2 by dim^2) stiffness tensor in material frame
  if (Dim == 1) {
    AKANTU_ERROR("Dimensions 1 not implemented: makes no sense to have "
                 "orthotropy for 1D");
  }

  Real Gamma;

  Gamma = 1 / (1 - _nu12 * _nu21 - _nu23 * _nu32 - _nu31 * _nu13 -
               2 * _nu21 * _nu32 * _nu13);

  // Lamé's first parameters
  _Cprime(0, 0) = _E1 * (1 - _nu23 * _nu32) * Gamma;
  _Cprime(1, 1) = _E2 * (1 - _nu13 * _nu31) * Gamma;
  _Cprime(1, 0) = _Cprime(0, 1) = _E1 * (_nu21 + _nu31 * _nu23) * Gamma;

  _Cprime(2, 2) = _E3 * (1 - _nu12 * _nu21) * Gamma;
  _Cprime(2, 0) = _Cprime(0, 2) = _E1 * (_nu31 + _nu21 * _nu32) * Gamma;
  _Cprime(2, 1) = _Cprime(1, 2) = _E2 * (_nu32 + _nu12 * _nu31) * Gamma;

  _Cprime(3, 3) = _G23;
  _Cprime(4, 4) = _G13;
  _Cprime(5, 5) = _G12;

  this->rotateCprime(_Cprime, _dir_vecs, _C);
  _C.eig(_eigC);
}

/* --------------------------------------------------------------------------
 */
template <>
void MaterialElasticOrthotropicHeterogeneous<2>::updateInternalParametersOnQuad(
    const Real & _E1, const Real & _E2, const Real & _E3, const Real & _nu12,
    const Real & _nu13, const Real & _nu23, const Real & _G12,
    const Real & /*_G13*/, const Real & /*_G23*/, Matrix<Real> & _Cprime,
    Matrix<Real> & _C, Vector<Real> & _eigC, Matrix<Real> & _dir_vecs) {

  Real _nu21 = _nu12 * _E2 / _E1;
  Real _nu31 = _nu13 * _E3 / _E1;
  Real _nu32 = _nu23 * _E3 / _E2;

  Real Gamma;

  if (this->plane_stress)
    Gamma = 1 / (1 - _nu12 * _nu21);
  else
    Gamma = 1 / (1 - _nu12 * _nu21 - _nu23 * _nu32 - _nu31 * _nu13 -
                 2 * _nu21 * _nu32 * _nu13);

  if (this->plane_stress) {
    _Cprime(0, 0) = _E1 * Gamma;
    _Cprime(1, 1) = _E2 * Gamma;
    _Cprime(1, 0) = _Cprime(0, 1) = _E1 * _nu21 * Gamma;
  } else {
    _Cprime(0, 0) = _E1 * (1 - _nu23 * _nu32) * Gamma;
    _Cprime(1, 1) = _E2 * (1 - _nu13 * _nu31) * Gamma;
    _Cprime(1, 0) = _Cprime(0, 1) = _E1 * (_nu21 + _nu31 * _nu23) * Gamma;
  }
  _Cprime(2, 2) = _G12;

  this->rotateCprime(_Cprime, _dir_vecs, _C);
  _C.eig(_eigC);
}

/* --------------------------------------------------------------------------
 */
template <UInt Dim>
void MaterialElasticOrthotropicHeterogeneous<
    Dim>::computePotentialEnergyByElement(ElementType type, UInt index,
                                          Vector<Real> & epot_on_quad_points) {

  Array<Real>::matrix_iterator gradu_it = this->gradu(type).begin(Dim, Dim);
  Array<Real>::matrix_iterator gradu_end = this->gradu(type).begin(Dim, Dim);
  Array<Real>::matrix_iterator stress_it = this->stress(type).begin(Dim, Dim);

  UInt nb_quadrature_points = this->fem.getNbIntegrationPoints(type);

  gradu_it += index * nb_quadrature_points;
  gradu_end += (index + 1) * nb_quadrature_points;
  stress_it += index * nb_quadrature_points;

  Real * epot_quad = epot_on_quad_points.storage();

  Matrix<Real> grad_u(Dim, Dim);

  for (; gradu_it != gradu_end; ++gradu_it, ++stress_it, ++epot_quad) {
    grad_u.copy(*gradu_it);

    this->computePotentialEnergyOnQuad(grad_u, *stress_it, *epot_quad);
  }
}

/* --------------------------------------------------------------------------
 */
INSTANTIATE_MATERIAL(elastic_orthotropic_heterogeneous,
                     MaterialElasticOrthotropicHeterogeneous);

} // namespace akantu
