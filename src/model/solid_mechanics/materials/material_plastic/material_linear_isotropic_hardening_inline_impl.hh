/**
 * @file   material_linear_isotropic_hardening_inline_impl.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Benjamin Paccaud <benjamin.paccaud@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Apr 07 2014
 * @date last modification: Thu Feb 20 2020
 *
 * @brief  Implementation of the inline functions of the material plasticity
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

#include "material_linear_isotropic_hardening.hh"

namespace akantu {
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/// Infinitesimal deformations
template <Int dim>
template <class Args, std::enable_if_t<not tuple::has_t<"F"_h, Args>::value> *>
inline void
MaterialLinearIsotropicHardening<dim>::computeStressOnQuad(Args && args) {

  Real delta_sigma_th = 0.;
  static_if(tuple::has_t<"sigma_th"_h, Args>())
      .then([&delta_sigma_th](auto && args) {
        delta_sigma_th = tuple::get<"previous_sigma_th"_h>(args) -
                         tuple::get<"sigma_th"_h>(args);
      })(std::forward<Args>(args));

  auto && grad_u = tuple::get<"grad_u"_h>(args);
  auto && previous_grad_u = tuple::get<"previous_grad_u"_h>(args);
  auto && grad_delta_u = grad_u - previous_grad_u;

  Matrix<Real, dim, dim> sigma_tr;

  // Compute trial stress, sigma_tr
  MaterialElastic<dim>::computeStressOnQuad(
      tuple::make_named_tuple(tuple::get<"grad_u"_h>() = grad_delta_u,
                              tuple::get<"sigma"_h>() = sigma_tr,
                              tuple::get<"sigma_th"_h>() = delta_sigma_th));

  auto && previous_sigma = tuple::get<"previous_sigma"_h>(args);
  sigma_tr += previous_sigma;

  // We need a full stress tensor, otherwise the VM stress is messed up
  Matrix<Real, 3, 3> sigma_tr_dev = Matrix<Real, 3, 3>::Zero();
  sigma_tr_dev.block<dim, dim>(0, 0) = sigma_tr;

  sigma_tr_dev -= Matrix<Real, 3, 3>::Identity() * sigma_tr.trace() / 3.0;

  // Compute effective deviatoric trial stress
  auto s = sigma_tr_dev.doubleDot(sigma_tr_dev);
  auto sigma_tr_dev_eff = std::sqrt(3. / 2. * s);

  auto && iso_hardening = tuple::get<"iso_hardening"_h>(args);
  auto initial_yielding =
      ((sigma_tr_dev_eff - iso_hardening - this->sigma_y) > 0);

  auto && previous_iso_hardening = tuple::get<"previous_iso_hardening"_h>(args);
  auto dp = (initial_yielding)
                ? (sigma_tr_dev_eff - this->sigma_y - previous_iso_hardening) /
                      (3 * this->mu + this->h)
                : 0;

  iso_hardening = previous_iso_hardening + this->h * dp;

  // Compute inelastic strain (ignore last components in 1-2D)
  Matrix<Real, dim, dim> delta_inelastic_strain =
      Matrix<Real, dim, dim>::Zero();
  if (std::abs(sigma_tr_dev_eff) > sigma_tr_dev.norm() * Math::getTolerance()) {
    delta_inelastic_strain =
        sigma_tr_dev.block<dim, dim>(0, 0) * 3. / 2. * dp / sigma_tr_dev_eff;
  }

  MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(tuple::append(
      args, tuple::get<"delta_inelastic_strain"_h>() = delta_inelastic_strain));
}

/* -------------------------------------------------------------------------- */
/// Finite deformations
template <Int dim>
template <class Args, std::enable_if_t<tuple::has_t<"F"_h, Args>::value> *>
inline void
MaterialLinearIsotropicHardening<dim>::computeStressOnQuad(Args && args) {
  // Finite plasticity
  Real dp = 0.0;
  Real d_dp = 0.0;
  UInt n = 0;

  Real delta_sigma_th = 0.;
  static_if(tuple::has_t<"sigma_th"_h, Args>())
      .then([&delta_sigma_th](auto && args) {
        delta_sigma_th = tuple::get<"previous_sigma_th"_h>(args) -
                         tuple::get<"sigma_th"_h>(args);
      })(std::forward<Args>(args));

  auto && grad_u = tuple::get<"grad_u"_h>(args);
  auto && previous_grad_u = tuple::get<"previous_grad_u"_h>(args);
  auto && grad_delta_u = grad_u - previous_grad_u;

  // Compute trial stress, sigma_tr
  Matrix<Real, dim, dim> sigma_tr;
  MaterialElastic<dim>::computeStressOnQuad(
      tuple::make_named_tuple(tuple::get<"grad_u"_h>() = grad_delta_u,
                              tuple::get<"sigma"_h>() = sigma_tr,
                              tuple::get<"sigma_th"_h>() = delta_sigma_th));

  auto && previous_sigma = tuple::get<"previous_sigma"_h>(args);
  sigma_tr += previous_sigma;

  // Compute deviatoric trial stress,  sigma_tr_dev
  auto && sigma_tr_dev =
      sigma_tr - Matrix<Real, dim, dim>::Identity() * sigma_tr.trace() / 3.0;

  // Compute effective deviatoric trial stress
  Real s = sigma_tr_dev.doubleDot(sigma_tr_dev);
  Real sigma_tr_dev_eff = std::sqrt(3. / 2. * s);

  auto && F = tuple::get<"F"_h>(args);

  // compute the cauchy stress to apply the Von-Mises criterion
  auto cauchy_stress = Material::StoCauchy<dim>(F, sigma_tr);

  auto && cauchy_stress_dev =
      cauchy_stress -
      Matrix<Real, dim, dim>::Identity() * cauchy_stress.trace() / 3.0;
  Real c = cauchy_stress_dev.doubleDot(cauchy_stress_dev);
  Real cauchy_stress_dev_eff = std::sqrt(3. / 2. * c);

  auto && iso_hardening = tuple::get<"iso_hardening"_h>(args);
  auto && previous_iso_hardening = tuple::get<"previous_iso_hardening"_h>(args);

  const Real iso_hardening_t = previous_iso_hardening;
  iso_hardening = iso_hardening_t;

  // F is written in terms of cauchy stress
  bool initial_yielding =
      ((cauchy_stress_dev_eff - iso_hardening - this->sigma_y) > 0);
  while (initial_yielding && std::abs(cauchy_stress_dev_eff - iso_hardening -
                                      this->sigma_y) > Math::getTolerance()) {

    d_dp = (cauchy_stress_dev_eff - 3. * this->mu * dp - iso_hardening -
            this->sigma_y) /
           (3. * this->mu + this->h);

    // r = r +  h * dp;
    dp = dp + d_dp;
    iso_hardening = iso_hardening_t + this->h * dp;

    ++n;
    /// TODO : explicit this criterion with an error message
    if ((d_dp < 1e-5) || (n > 50)) {
      AKANTU_DEBUG_INFO("convergence of increment of plastic strain. d_dp:"
                        << d_dp << "\tNumber of iteration:" << n);
      break;
    }
  }

  // Update internal variable
  Matrix<Real, dim, dim> delta_inelastic_strain;
  if (std::abs(sigma_tr_dev_eff) > sigma_tr_dev.norm() * Math::getTolerance()) {

    // /// compute the direction of the plastic strain as \frac{\partial
    // F}{\partial S} = \frac{3}{2J\sigma_{effective}}} Ft \sigma_{dev} F
    Matrix<Real, dim, dim> cauchy_dev_F;
    cauchy_dev_F = F * cauchy_stress_dev;
    Real J = F.determinant();
    Real constant = not Math::are_float_equal(J, 0.) ? 1. / J : 0;
    constant *= 3. * dp / (2. * cauchy_stress_dev_eff);
    delta_inelastic_strain = F.transpose() * cauchy_dev_F * constant;

    // Direction given by the piola kirchhoff deviatoric tensor \frac{\partial
    // F}{\partial S} = \frac{3}{2\sigma_{effective}}}S_{dev}
    // delta_inelastic_strain.copy(sigma_tr_dev);
    // delta_inelastic_strain *= 3./2. * dp / sigma_tr_dev_eff;
  }

  MaterialPlastic<dim>::computeStressAndInelasticStrainOnQuad(tuple::append(
      args, tuple::get<"delta_inelastic_strain"_h>() = delta_inelastic_strain));
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void MaterialLinearIsotropicHardening<dim>::computeTangentModuliOnQuad(
    Args && args) const {
  // Initial tangent
  MaterialElastic<dim>::computeTangentModuliOnQuad(args);
}
} // namespace akantu
