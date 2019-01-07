/**
 * @file   material_damage_iterative_viscoelastic.cc
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date   Tue Nov 20 2016
 *
 * @brief  Inline implementation of material iterative stiffness reduction viscoelastic
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
#include "material_damage_iterative_viscoelastic.hh"
#include "communicator.hh"
#include "solid_mechanics_model_RVE.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, Real & dam) {

  Real dt = this->model.getTimeStep();
  Real E_ef = this->Einf * (1 - dam);

  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Real lambda = (1 - dam) * this->Eta(k) / this->Ev(k) / (1 - dam);
    Real exp_dt_lambda = std::exp(-dt / lambda);

    if (exp_dt_lambda == 1) {
      E_ef += (1 - dam) * this->Ev(k);
    } else {
      E_ef += std::min((1 - exp_dt_lambda) * (1 - dam) * this->Ev(k) * lambda / dt, (1 - dam) * this->Ev(k)); ;
    }
  }

  tangent.copy(this->C);
  tangent *= E_ef;

}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::updateIntVarOnQuad(
    Matrix<Real> grad_u, Matrix<Real> previous_grad_u,
    Tensor3<Real> & sigma_v, Tensor3<Real> & epsilon_v, Real dam) {

  Matrix<Real> grad_delta_u(grad_u);
  grad_delta_u -= previous_grad_u;

  Real dt = this->model.getTimeStep();
  Vector<Real> voigt_delta_strain(voigt_h::size);
  for (UInt I = 0; I < voigt_h::size; ++I) {
      Real voigt_factor = voigt_h::factors[I];
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      voigt_delta_strain(I) =
          voigt_factor * (grad_delta_u(i, j) + grad_delta_u(j, i)) / 2.;
  }


  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Real lambda = (1 - dam) * this->Eta(k) / this->Ev(k) / (1 - dam);
    Real exp_dt_lambda = exp(-dt / lambda);
    Real E_ef_v;

    if (exp_dt_lambda == 1) {
      E_ef_v = (1 - dam) * this->Ev(k);
    } else {
      E_ef_v = std::min((1 - exp_dt_lambda) * (1 - dam) * this->Ev(k) * lambda / dt, (1 - dam) * this->Ev(k));
    }

    Vector<Real> voigt_sigma_v(voigt_h::size);
    Vector<Real> voigt_epsilon_v(voigt_h::size);

    for (UInt I = 0; I < voigt_h::size; ++I) {
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      voigt_sigma_v(I) = sigma_v(i, j, k);
    }

    voigt_sigma_v = exp_dt_lambda * voigt_sigma_v + E_ef_v * this->C *
        voigt_delta_strain;
    voigt_epsilon_v = 1 / this->Ev(k) / (1 - dam) *
        this->S * voigt_sigma_v;

    for (UInt I = 0; I < voigt_h::size; ++I) {
      Real voigt_factor = voigt_h::factors[I];
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      sigma_v(i, j, k) = sigma_v(j, i, k) = voigt_sigma_v(I);
      epsilon_v(i, j, k) = epsilon_v(j, i, k) =
          voigt_epsilon_v(I) / voigt_factor;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::updateDissipatedEnergyDamageOnQuad( Matrix<Real> grad_u, Matrix<Real> epsilon_p, Tensor3<Real> sigma_v, Tensor3<Real> epsilon_v, Tensor3<Real> sigma_v_pr, Tensor3<Real> epsilon_v_pr, bool damaged, Real dam, Real dam_pr, Real & epot, Real & ints, Real & edd) {

  Matrix<Real> delta_grad_u(grad_u);
  delta_grad_u -= epsilon_p;

  Matrix<Real> sigma_h(spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_p(spatial_dimension, spatial_dimension);

  Real trace_h = grad_u.trace(); // trace = (\nabla u)_{kk}
  Real trace_p = epsilon_p.trace();

  // \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  // u_{ij} + \nabla u_{ji})
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      sigma_h(i, j) = (1 - dam) * ((i == j) * this->lambda * trace_h +
                                   this->mu * (grad_u(i, j) +
                                               grad_u(j, i)));
      sigma_p(i, j) = (1 - dam_pr) * ((i == j) * this->lambda * trace_p +
                                      this->mu * (epsilon_p(i, j) +
                                                  epsilon_p(j, i)));
    }
  }

  sigma_h += sigma_p;

  Real dint = .5 * sigma_h.doubleDot(delta_grad_u);

  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Matrix<Real> stress_v_h(sigma_v(k));
    Matrix<Real> stress_v_pr(sigma_v_pr(k));
    stress_v_h += stress_v_pr;
    Matrix<Real> delta_strain_v(epsilon_v(k));
    Matrix<Real> strain_v_pr(epsilon_v_pr(k));
    delta_strain_v -= strain_v_pr;

    dint += .5 * stress_v_h.doubleDot(delta_strain_v);
  }

  /// update elastic part of the mechanical work for all the elements
  ints += dint;

  /// update dissipated energy only from freshly damaged elements
  if (damaged) {
    edd = ints - epot;
    if (std::abs(edd) < std::numeric_limits<Real>::epsilon())
      edd = 0;
  }
}


/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::computeStressOnQuad(
    Matrix<Real> grad_u, Matrix<Real> previous_grad_u,
    Matrix<Real> & sigma, Tensor3<Real> sigma_v, Real sigma_th, Real damage) {

  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  Vector<Real> voigt_current_strain(voigt_h::size);
  Vector<Real> voigt_previous_strain(voigt_h::size);
  Vector<Real> voigt_stress(voigt_h::size);
  Vector<Real> voigt_sigma_v(voigt_h::size);

  for (UInt I = 0; I < voigt_h::size; ++I) {
    Real voigt_factor = voigt_h::factors[I];
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    voigt_current_strain(I) =
        voigt_factor * (grad_u(i, j) + grad_u(j, i)) / 2.;
    voigt_previous_strain(I) =
        voigt_factor * (previous_grad_u(i, j) + previous_grad_u(j, i)) / 2.;
  }

  voigt_stress = (1 - damage) * this->Einf * this->C * voigt_current_strain;

  Real dt = this->model.getTimeStep();

  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Real lambda = (1 - damage) * this->Eta(k) / this->Ev(k)/ (1 - damage);
    Real exp_dt_lambda = exp(-dt / lambda);
    Real E_additional;

    if (exp_dt_lambda == 1) {
      E_additional = (1 - damage) * this->Ev(k);
    } else {
      E_additional = std::min((1 - exp_dt_lambda) * (1 - damage) * this->Ev(k) * lambda / dt, (1 - damage) * this->Ev(k));
    }

    for (UInt I = 0; I < voigt_h::size; ++I) {
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      voigt_sigma_v(I) = sigma_v(i, j, k);
    }

    voigt_stress += E_additional * this->C * (voigt_current_strain - voigt_previous_strain) +
        exp_dt_lambda * voigt_sigma_v;
  }

  for (UInt I = 0; I < voigt_h::size; ++I) {
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    sigma(i, j) = sigma(j, i) =
        voigt_stress(I) + (i == j) * sigma_th;
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::computePotentialEnergyOnQuad(
    Matrix<Real> grad_u, Real & epot,
    Tensor3<Real> sigma_v, Tensor3<Real> epsilon_v, Real dam) {

  Real trace = grad_u.trace(); // trace = (\nabla u)_{kk}

  Matrix<Real> sigma(spatial_dimension, spatial_dimension);
  // \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  // u_{ij} + \nabla u_{ji})
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      sigma(i, j) = (1 - dam) * ((i == j) * this->lambda * trace +
                                 this->mu * (grad_u(i, j) +
                                             grad_u(j, i)));
    }
  }

  epot = 0.5 * sigma.doubleDot(grad_u);

  for (UInt k = 0; k < this->Eta.size(); ++k) {
    Matrix<Real> stress_v = sigma_v(k);
    Matrix<Real> strain_v = epsilon_v(k);
    epot += 0.5 * stress_v.doubleDot(strain_v);
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
