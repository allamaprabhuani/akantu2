/**
 * @file   weight_function.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Mar 23 15:55:58 2012
 *
 * @brief  implementation of the weight function classes
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

/* -------------------------------------------------------------------------- */
/* Stress based weight function                                               */
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
StressBasedWeightFunction<spatial_dimension>::StressBasedWeightFunction(const Material & material) :
  BaseWeightFunction<spatial_dimension>(material),
  ft(0.),
  stress_diag("stress_diag", material.getID()), selected_stress_diag(NULL),
  stress_base("stress_base", material.getID()), selected_stress_base(NULL),
  characteristic_size("lc", material.getID()),  selected_characteristic_size(NULL)
{
  material.initInternalVector(stress_diag, spatial_dimension);
  material.initInternalVector(stress_base, spatial_dimension * spatial_dimension);
  material.initInternalVector(characteristic_size, 1);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void StressBasedWeightFunction<spatial_dimension>::init() {
  this->material.resizeInternalVector(stress_diag);
  this->material.resizeInternalVector(stress_base);
  this->material.resizeInternalVector(characteristic_size);

  const Mesh & mesh = this->material.getModel().getFEM().getMesh();
  for (UInt g = _not_ghost; g < _casper; ++g) {
    GhostType gt = GhostType(g);
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt);
    for(; it != last_type; ++it) {
      UInt nb_quadrature_points =
	this->material.getModel().getFEM().getNbQuadraturePoints(*it, gt);
      const Vector<UInt> & element_filter = this->material.getElementFilter(*it, gt);
      UInt nb_element = element_filter.getSize();

      Vector<Real> ones(nb_element*nb_quadrature_points, 1, 1.);
      Vector<Real> & lc = characteristic_size(*it, gt);
      this->material.getModel().getFEM().integrateOnQuadraturePoints(ones,
								     lc,
								     1,
								     *it,
								     gt,
								     &element_filter);

      for (UInt q = 0;  q < nb_quadrature_points * nb_element; q++) {
	lc(q) = pow(lc(q), 1./ Real(spatial_dimension));
      }
    }
  }
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void StressBasedWeightFunction<spatial_dimension>::updatePrincipalStress(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = this->material.getModel().getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    Vector<Real>::const_iterator<types::Matrix> sigma =
      this->material.getStress(*it, ghost_type).begin(spatial_dimension, spatial_dimension);
    Vector<Real>::iterator<types::RVector> eigenvalues =
      stress_diag(*it, ghost_type).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> eigenvalues_end =
      stress_diag(*it, ghost_type).end(spatial_dimension);
    Vector<Real>::iterator<types::Matrix> eigenvector =
      stress_base(*it, ghost_type).begin(spatial_dimension, spatial_dimension);

#ifndef __trick__
    Vector<Real>::iterator<Real> cl = characteristic_size(*it, ghost_type).begin();
#endif
    UInt q = 0;
    for(;eigenvalues != eigenvalues_end; ++sigma, ++eigenvalues, ++eigenvector, ++cl, ++q) {
      //      if (q == 4871) std::cout << "Sigma " << *sigma;
      sigma->eig(*eigenvalues, *eigenvector);
      //      if (q == 17774) std::cout << "Before " << *eigenvalues;
      *eigenvalues /= ft;
#ifndef __trick__
      // specify a lower bound for principal stress based on the size of the element
      for (UInt i = 0; i < spatial_dimension; ++i) {
        (*eigenvalues)(i) = std::max(*cl / this->R, (*eigenvalues)(i));
      }
#endif
      //      if (q == 17774) std::cout << "After " << *eigenvalues;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void StressBasedWeightFunction<spatial_dimension>::selectType(ElementType type1,
								     GhostType ghost_type1,
								     ElementType type2,
								     GhostType ghost_type2) {
  selected_stress_diag = &stress_diag(type2, ghost_type2);
  selected_stress_base = &stress_base(type2, ghost_type2);

  selected_characteristic_size = &characteristic_size(type1, ghost_type1);
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real StressBasedWeightFunction<spatial_dimension>::operator()(Real r,
								     QuadraturePoint & q1,
								     QuadraturePoint & q2) {
  Real zero = std::numeric_limits<Real>::epsilon();

  if(r < zero) return 1.; // means x and s are the same points

  const types::RVector & x = q1.getPosition();
  const types::RVector & s = q2.getPosition();

  types::RVector eigs =
    selected_stress_diag->begin(spatial_dimension)[q2.global_num];

  types::Matrix eigenvects =
    selected_stress_base->begin(spatial_dimension, spatial_dimension)[q2.global_num];

  Real min_rho_lc = selected_characteristic_size->begin()[q1.global_num];

  types::RVector x_s(spatial_dimension);
  x_s  = x;
  x_s -= s;

  Real rho_2 = 0.;

  switch(spatial_dimension) {
  case 1: {
    rho_2 = eigs[0];
    break;
  }
  case 2:{
    types::RVector u1(eigenvects.storage(), spatial_dimension);
    Real cos_t = x_s.dot(u1) / (x_s.norm() * u1.norm());

    Real cos_t_2;
    Real sin_t_2;

    Real sigma1_2 = eigs[0]*eigs[0];
    Real sigma2_2 = eigs[1]*eigs[1];

#ifdef __trick__
    if(std::abs(cos_t) < zero) {
      cos_t_2 = 0;
      sin_t_2 = 1;
    } else {
#endif
      cos_t_2 = cos_t * cos_t;
      sin_t_2 = (1 - cos_t_2);
#ifdef __trick__
    }

    Real rhop1 = std::max(0., cos_t_2 / sigma1_2);
    Real rhop2 = std::max(0., sin_t_2 / sigma2_2);
#else
    Real rhop1 = cos_t_2 / sigma1_2;
    Real rhop2 = sin_t_2 / sigma2_2;
#endif

    rho_2 = 1./ (rhop1 + rhop2);

    break;
  }
  case 3: {
    types::RVector u1(eigenvects.storage(), spatial_dimension);
    types::RVector u3(eigenvects.storage() + 2*spatial_dimension, spatial_dimension);

    types::RVector tmp(spatial_dimension);
    tmp.crossProduct(x_s, u3);

    types::RVector u3_C_x_s_C_u3(spatial_dimension);
    u3_C_x_s_C_u3.crossProduct(u3, tmp);

    Real norm_u3_C_x_s_C_u3 = u3_C_x_s_C_u3.norm();
    Real cos_t = 0.;
    if(std::abs(norm_u3_C_x_s_C_u3) > zero) {
      Real inv_norm_u3_C_x_s_C_u3 = 1. / norm_u3_C_x_s_C_u3;
      cos_t = u1.dot(u3_C_x_s_C_u3) * inv_norm_u3_C_x_s_C_u3;
    }

    Real cos_p = u3.dot(x_s) / r;
    //    std::cout << " - cos_t: " << std::setw(13) << cos_t;
    //    std::cout << " - cos_p: " << std::setw(13) << cos_p;

    Real cos_t_2;
    Real sin_t_2;
    Real cos_p_2;
    Real sin_p_2;

    Real sigma1_2 = eigs[0]*eigs[0];
    Real sigma2_2 = eigs[1]*eigs[1];
    Real sigma3_2 = eigs[2]*eigs[2];

#ifdef __trick__
    if(std::abs(cos_t) < zero) {
      cos_t_2 = 0;
      sin_t_2 = 1;
    } else {
#endif
      cos_t_2 = cos_t * cos_t;
      sin_t_2 = (1 - cos_t_2);
#ifdef __trick__
    }

    if(std::abs(cos_p) < zero) {
      cos_p_2 = 0;
      sin_p_2 = 1;
    } else {
#endif
      cos_p_2 = cos_p * cos_p;
      sin_p_2 = (1 - cos_p_2);
#ifdef __trick__
    }

    Real rhop1 = std::max(0., sin_p_2 * cos_t_2 / sigma1_2);
    Real rhop2 = std::max(0., sin_p_2 * sin_t_2 / sigma2_2);
    Real rhop3 = std::max(0., cos_p_2 / sigma3_2);
#else
    Real rhop1 = sin_p_2 * cos_t_2 / sigma1_2;
    Real rhop2 = sin_p_2 * sin_t_2 / sigma2_2;
    Real rhop3 = cos_p_2 / sigma3_2;
#endif

    rho_2 = 1./ (rhop1 + rhop2 + rhop3);
  }
  }

  Real rho_lc_2 = std::max(this->R2 * rho_2, min_rho_lc*min_rho_lc);

  Real w = std::max(0., 1. - r*r / rho_lc_2);
  w = w*w;

  //    std::cout << "(" << q1 << "," << q2 << ") " << w << std::endl;
  //  Real w = exp(- 2*2*r*r / rho_lc_2);

  return w;
}
