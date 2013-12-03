/**
 * @file   material_elastic_linear_anisotropic.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   Tue May 08 13:01:18 2012
 *
 * @brief  Anisotropic elastic material
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
#include "material_elastic_linear_anisotropic.hh"
#include "solid_mechanics_model.hh"
#include <algorithm>
#include <sstream>


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialElasticLinearAnisotropic<spatial_dimension>::
MaterialElasticLinearAnisotropic(SolidMechanicsModel & model,
                                 const ID & id,
                                 bool symmetric)  :
  Material(model, id),
  rot_mat(spatial_dimension, spatial_dimension),
  Cprime(spatial_dimension*spatial_dimension,
         spatial_dimension*spatial_dimension),
  C(this->voigt_h.size, this->voigt_h.size),
  eigC(this->voigt_h.size),
  symmetric(symmetric),
  alpha(0){
  AKANTU_DEBUG_IN();

  this->dir_vecs.push_back(new Vector<Real>(spatial_dimension));
  (*this->dir_vecs.back())[0] = 1.;
  this->registerParam("n1", *(this->dir_vecs.back()), _pat_parsmod,
                      "Direction of main material axis");
  this->dir_vecs.push_back(new Vector<Real>(spatial_dimension));
  (*this->dir_vecs.back())[1] = 1.;
  this->registerParam("n2", *(this->dir_vecs.back()), _pat_parsmod,
                      "Direction of secondary material axis");

  if (spatial_dimension > 2) {
    this->dir_vecs.push_back(new Vector<Real>(spatial_dimension));
    (*this->dir_vecs.back())[2] = 1.;
    this->registerParam("n3", *(this->dir_vecs.back()), _pat_parsmod,
                        "Direction of tertiary material axis");
  }

  for (UInt i = 0 ;  i < this->voigt_h.size ; ++i) {
    UInt start = 0;
    if (this->symmetric) {
      start = i;
    }
    for (UInt j = start ;  j < this->voigt_h.size ; ++j) {
      std::stringstream param("C");
      param << "C" << i+1 << j+1;
      this->registerParam(param.str() , this->Cprime(i,j), 0., _pat_parsmod,
                          "Coefficient " + param.str());
    }
  }

  this->registerParam("alpha", this->alpha, _pat_parsmod,
                      "Proportion of viscous stress");

  AKANTU_DEBUG_OUT();
  }

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialElasticLinearAnisotropic<spatial_dimension>::~MaterialElasticLinearAnisotropic() {
  for (UInt i = 0 ;  i < spatial_dimension ; ++i) {
    delete this->dir_vecs[i];
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticLinearAnisotropic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticLinearAnisotropic<spatial_dimension>::updateInternalParameters() {
  Material::updateInternalParameters();
  if (this->symmetric) {
    for (UInt i = 0 ;  i < this->voigt_h.size ; ++i) {
      for (UInt j = i+1 ;  j < this->voigt_h.size ; ++j) {
        this->Cprime(j, i) = this->Cprime(i, j);
      }
    }
  }
  this->rotateCprime();
  this->C.eig(this->eigC);
}

/* -------------------------------------------------------------------------- */
template<UInt Dim>
void MaterialElasticLinearAnisotropic<Dim>::rotateCprime() {
  // start by filling the empty parts fo Cprime
  UInt diff = Dim*Dim - this->voigt_h.size;
  for (UInt i = this->voigt_h.size ;  i < Dim*Dim ; ++i) {
    for (UInt j = 0 ;  j < Dim*Dim ; ++j) {
      this->Cprime(i, j) = this->Cprime(i-diff, j);
    }
  }
  for (UInt i = 0 ;  i < Dim*Dim ; ++i) {
    for (UInt j = this->voigt_h.size ;  j < Dim*Dim ; ++j) {
      this->Cprime(i, j) = this->Cprime(i, j-diff);
    }
  }
  // construction of rotator tensor
  // normalise rotation matrix
  for (UInt j = 0 ;  j < Dim ; ++j) {
    this->rot_mat(j) = *this->dir_vecs[j];
    this->rot_mat(j).normalize();
  }

  // make sure the vectors form a right-handed base
  Vector<Real> test_axis(3);
  Vector<Real> v1(3),v2(3),v3(3);

  if (Dim == 2){
    for (UInt i = 0; i < Dim; ++i) {
      v1[i] = this->rot_mat(0)[i];
      v2[i] = this->rot_mat(1)[i];
      v3[i] = 0.;
    }
    v3[2] = 1.;
    v1[2] = 0.;
    v2[2] = 0.;
  }
  else if (Dim == 3){
    v1 = this->rot_mat(0);
    v2 = this->rot_mat(1);
    v3 = this->rot_mat(2);
  }
  
  test_axis.crossProduct(v1,v2);
  test_axis -= v3;
  if (test_axis.norm() > 8*std::numeric_limits<Real>::epsilon()) {
    AKANTU_DEBUG_ERROR("The axis vectors do not form a right-handed coordinate "
                       << "system. I. e., ||n1 x n2 - n3|| should be zero, but "
                       << "it is " << test_axis.norm() << ".");
  }

  // create the rotator and the reverse rotator
  Matrix<Real> rotator(Dim * Dim, Dim * Dim);
  Matrix<Real> revrotor(Dim * Dim, Dim * Dim);
  for (UInt i = 0 ;  i < Dim ; ++i) {
    for (UInt j = 0 ;  j < Dim ; ++j) {
      for (UInt k = 0 ;  k < Dim ; ++k) {
        for (UInt l = 0 ;  l < Dim ; ++l) {
          UInt I = this->voigt_h.mat[i][j];
          UInt J = this->voigt_h.mat[k][l];
          rotator (I, J) = this->rot_mat(k, i) * this->rot_mat(l, j);
          revrotor(I, J) = this->rot_mat(i, k) * this->rot_mat(j, l);
        }
      }
    }
  }

  // create the full rotated matrix
  Matrix<Real> Cfull(Dim*Dim, Dim*Dim);
  Cfull = rotator*Cprime*revrotor;

  for (UInt i = 0 ;  i < this->voigt_h.size ; ++i) {
    for (UInt j = 0 ;  j < this->voigt_h.size ; ++j) {
      this->C(i, j) = Cfull(i, j);
    }
  }

}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticLinearAnisotropic<spatial_dimension>::
computeStress(ElementType el_type, GhostType ghost_type) {
  if (this->alpha == 0.) {
    this->computeStressWorker<false>(el_type, ghost_type);
  } else {
    this->computeStressWorker<true>(el_type, ghost_type);
  }
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
template<bool viscous>
void MaterialElasticLinearAnisotropic<spatial_dimension>::
computeStressWorker(ElementType el_type, GhostType ghost_type) {

  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  AKANTU_DEBUG_IN();


  Array<Real>::iterator< Matrix<Real> > strain_it =
    this->strain(el_type, ghost_type).begin(spatial_dimension,
                                            spatial_dimension);
  Array<Real>::iterator< Matrix<Real> > strain_end =
    this->strain(el_type, ghost_type).end(spatial_dimension,
                                          spatial_dimension);

  Array<Real> strain_rate(0, spatial_dimension * spatial_dimension,
                          "strain_rate");
  Array<Real>::iterator< Matrix<Real> > strain_rate_it;
  if (viscous) {
    const Array<UInt> & elem_filter = this->element_filter(el_type, ghost_type);
    Array<Real> & velocity = this->model->getVelocity();
    this->model->getFEM().gradientOnQuadraturePoints(velocity, strain_rate,
                                                     spatial_dimension, el_type,
                                                     ghost_type,
                                                     elem_filter);
    strain_rate_it = strain_rate.begin(spatial_dimension, spatial_dimension);
  }
  this->stress(el_type, ghost_type).resize(this->strain(el_type,
                                                        ghost_type).getSize());

  UInt nb_quad_pts = strain_end - strain_it;

  // create array for strains and stresses of all dof of all gauss points
  // for efficient computation of stress
  Matrix<Real> voigt_strains(this->voigt_h.size, nb_quad_pts);
  Matrix<Real> voigt_stresses(this->voigt_h.size, nb_quad_pts);
  Matrix<Real> * grad_u_dot;
  // copy strains
  for (UInt q = 0; strain_it != strain_end ; ++strain_it, ++q) {
    Matrix<Real> & grad_u = *strain_it;
    if (viscous) {
      grad_u_dot = &(*strain_rate_it);
    }
    for (UInt i = 0 ; i < spatial_dimension ; ++i) {
      for (UInt j = i ; j < spatial_dimension ; ++j) {
        UInt I = this->voigt_h.mat[i][j];
        voigt_strains(I, q) = 0.5 * this->voigt_h.factors[I]*(grad_u(j,i)
                                                              + grad_u(i,j));
        if (viscous) {
          voigt_strains(I, q) += (0.5 * this->alpha * this->voigt_h.factors[I] *
                                  (grad_u_dot->operator()(j, i)
                                   + grad_u_dot->operator()(i, j)));
        }
      }
    }
    if (viscous) {
      ++strain_rate_it;
    }
  }

  // compute stresses
  voigt_stresses = this->C * voigt_strains;

  // copy stresses back
  Array<Real>::iterator< Matrix<Real> > stress_it =
    this->stress(el_type, ghost_type).begin(spatial_dimension,
                                            spatial_dimension);
  Array<Real>::iterator< Matrix<Real> > stress_end =
    this->stress(el_type, ghost_type).end(spatial_dimension,
                                          spatial_dimension);
  for (UInt q = 0 ; stress_it != stress_end; ++stress_it, ++q) {
    Matrix<Real> & stress = *stress_it;
    for (UInt i = 0 ;  i < spatial_dimension ; ++i) {
      stress(i,i) = voigt_stresses(this->voigt_h.mat[i][i], q);
      for (UInt j = i+1 ;  j < spatial_dimension ; ++j) {
        UInt I = this->voigt_h.mat[i][j];
        stress(j, i) =
          stress(i, j) = voigt_stresses(I, q);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticLinearAnisotropic<spatial_dimension>::computeTangentModuli(const ElementType & el_type,
                                                                               Array<Real> & tangent_matrix,
                                                                               GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  tangent.copy(this->C);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<UInt dim>
Real MaterialElasticLinearAnisotropic<dim>::getStableTimeStep(Real h,
                                                              __attribute__ ((unused)) const Element & element) {
  return h/this->getCelerity();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElasticLinearAnisotropic<spatial_dimension>::getCelerity() const {
  return std::sqrt( this->eigC(0) / rho);
}


/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialElasticLinearAnisotropic);


__END_AKANTU__
