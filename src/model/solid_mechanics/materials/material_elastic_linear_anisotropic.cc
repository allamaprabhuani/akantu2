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
MaterialElasticLinearAnisotropic<spatial_dimension>::MaterialElasticLinearAnisotropic(
          SolidMechanicsModel & model,
          const ID & id,
          bool symmetric)  :
  Material(model, id),
  C(this->getTangentStiffnessVoigtSize(spatial_dimension),
    this->getTangentStiffnessVoigtSize(spatial_dimension)),
  symmetric(symmetric){
  AKANTU_DEBUG_IN();
  for (UInt i = 0 ;  i < this->voigt_h.size ; ++i) {
    UInt start = 0;
    if (this->symmetric) {
      start = i;
    }
    for (UInt j = start ;  j < this->voigt_h.size ; ++j) {
      std::stringstream param("C");
      param << i+1 << j+1;
      this->registerParam(param.str() , this->C(i,j), 0., _pat_parsmod,
                          "Coefficient " + param.str());
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialElasticLinearAnisotropic<spatial_dimension>::~MaterialElasticLinearAnisotropic() {
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
        this->C(j, i) = this->C(i, j);
      }
    }
  }
  Matrix <Real> trash_eigen_vectors;
  C.eig(this->eigC, trash_eigen_vectors);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticLinearAnisotropic<spatial_dimension>::computeStress(ElementType el_type,
								  GhostType ghost_type) {

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
  this->stress(el_type, ghost_type).resize(this->strain(el_type,
                                                        ghost_type).getSize());

  UInt nb_quad_pts = strain_end - strain_it;

  // create array for strains and stresses of all dof of all gauss points
  // for efficient computation of stress
  Matrix<Real> voigt_strains(this->voigt_h.size, nb_quad_pts);
  Matrix<Real> voigt_stresses(this->voigt_h.size, nb_quad_pts);

  // copy strains
  for (UInt q = 0; strain_it != strain_end ; ++strain_it, ++q) {
    Matrix<Real> & grad_u = *strain_it;
    for (UInt j = 0 ; j < spatial_dimension ; ++j) {
      for (UInt i = j ; i < spatial_dimension ; ++i) {
        UInt I = this->voigt_h.mat[i][j];
        voigt_strains(I, q) = this->voigt_h.factors[I]*(grad_u(j,i) + grad_u(i,j));
      }
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
    for (UInt j = 0 ;  j < spatial_dimension ; ++j) {
      stress(j,j) = voigt_stresses(this->voigt_h.mat[j][j], q);
      for (UInt i = j+1 ;  i < spatial_dimension ; ++i) {
        stress(j, i) =
          stress(i, j) = voigt_stresses(this->voigt_h.mat[i][j], q);
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
