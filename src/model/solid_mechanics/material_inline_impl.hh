/**
 * @file   material_inline_impl.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Implementation of the inline functions of the class material
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt Material::getTangentStiffnessVoigtSize(UInt dim) {
  return (dim * (dim - 1) / 2 + dim);
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getCauchyStressMatrixSize(UInt dim) {
  return (dim * dim);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::gradUToF(const Matrix<Real> & grad_u, Matrix<Real> & F) {
  AKANTU_DEBUG_ASSERT(F.size() >= grad_u.size() && grad_u.size() == dim * dim,
                      "The dimension of the tensor F should be greater or "
                      "equal to the dimension of the tensor grad_u.");
  F.eye();

  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j) {
      F(i, j) += grad_u(i, j);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline decltype(auto) Material::gradUToF(const Matrix<Real> & grad_u) {
  Matrix<Real> F(dim, dim);
  gradUToF<dim>(grad_u, F);
  return F;
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::StoCauchy(const Matrix<Real> & F, const Matrix<Real> & S,
                                Matrix<Real> & sigma, const Real & C33) const {
  Real J = F.det() * sqrt(C33);

  Matrix<Real> F_S(dim, dim);
  F_S = F * S;
  Real constant = J ? 1. / J : 0;
  sigma.mul<false, true>(F_S, F, constant);
}

/* -------------------------------------------------------------------------- */
inline void Material::rightCauchy(const Matrix<Real> & F, Matrix<Real> & C) {
  C.mul<true, false>(F, F);
}

/* -------------------------------------------------------------------------- */
inline void Material::leftCauchy(const Matrix<Real> & F, Matrix<Real> & B) {
  B.mul<false, true>(F, F);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::gradUToEpsilon(const Matrix<Real> & grad_u,
                                     Matrix<Real> & epsilon) {
  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j) {
      epsilon(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline decltype(auto) Material::gradUToEpsilon(const Matrix<Real> & grad_u) {
  Matrix<Real> epsilon(dim, dim);
  Material::template gradUToEpsilon<dim>(grad_u, epsilon);
  return epsilon;
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::gradUToE(const Matrix<Real> & grad_u, Matrix<Real> & E) {
  E.mul<true, false>(grad_u, grad_u, .5);

  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j) {
      E(i, j) += 0.5 * (grad_u(i, j) + grad_u(j, i));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline decltype(auto) Material::gradUToE(const Matrix<Real> & grad_u) {
  Matrix<Real> E(dim, dim);
  gradUToE<dim>(grad_u, E);
  return E;
}

/* -------------------------------------------------------------------------- */
inline Real Material::stressToVonMises(const Matrix<Real> & stress) {
  // compute deviatoric stress
  UInt dim = stress.cols();
  Matrix<Real> deviatoric_stress =
      Matrix<Real>::eye(dim, -1. * stress.trace() / 3.);

  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j) {
      deviatoric_stress(i, j) += stress(i, j);
    }
  }

  // return Von Mises stress
  return std::sqrt(3. * deviatoric_stress.doubleDot(deviatoric_stress) / 2.);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::setCauchyStressMatrix(const Matrix<Real> & S_t,
                                            Matrix<Real> & sigma) {
  AKANTU_DEBUG_IN();

  sigma.zero();

  /// see Finite ekement formulations for large deformation dynamic analysis,
  /// Bathe et al. IJNME vol 9, 1975, page 364 ^t \f$\tau\f$
  for (UInt i = 0; i < dim; ++i) {
    for (UInt m = 0; m < dim; ++m) {
      for (UInt n = 0; n < dim; ++n) {
        sigma(i * dim + m, i * dim + n) = S_t(m, n);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getNbData(const Array<Element> & elements,
                                const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    return (this->isFiniteDeformation() ? 3 : 1) * spatial_dimension *
           spatial_dimension * sizeof(Real) *
           this->getModel().getNbIntegrationPoints(elements);
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void Material::packData(CommunicationBuffer & buffer,
                               const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      packInternalFieldHelper(piola_kirchhoff_2, buffer, elements);
      packInternalFieldHelper(gradu, buffer, elements);
    }
    packInternalFieldHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::unpackData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      unpackInternalFieldHelper(piola_kirchhoff_2, buffer, elements);
      unpackInternalFieldHelper(gradu, buffer, elements);
    }
    unpackInternalFieldHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Material::inflateInternal(const std::string & field_id,
                               const ElementTypeMapArray<T> & field,
                               GhostType ghost_type, ElementKind element_kind) {
  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in material "
                                                   << this->name);
  }

  InternalField<T> & internal_field = this->template getInternal<T>(field_id);
  const FEEngine & fe_engine = internal_field.getFEEngine();

  for (auto && type : field.elementTypes(ghost_type)) {
    if (not internal_field.exists(type, ghost_type)) {
      continue;
    }
    const auto & filter = internal_field.getFilter(type, ghost_type);

    const auto & src_array = field(type, ghost_type);
    auto & dest_array = internal_field(type, ghost_type);

    auto nb_quad_per_elem = fe_engine.getNbIntegrationPoints(type);
    auto nb_component = src_array.getNbComponent();

    AKANTU_DEBUG_ASSERT(
        field.size() == fe_engine.getMesh().getNbElement(type, ghost_type) *
                            nb_quad_per_elem,
        "The ElementTypeMapArray to inflate is not of the proper size");
    AKANTU_DEBUG_ASSERT(
        dest_array.getNbComponent() == nb_component,
        "The ElementTypeMapArray has not the proper number of components");

    auto src =
        make_view(field(type, ghost_type), nb_component, nb_quad_per_elem)
            .begin();
    for (auto && data :
         zip(filter, make_view(dest_array, nb_component, nb_quad_per_elem))) {
      std::get<1>(data) = src[std::get<0>(data)];
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const InternalField<T> &
Material::getInternal([[gnu::unused]] const ID & int_id) const {
  AKANTU_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline InternalField<T> &
Material::getInternal([[gnu::unused]] const ID & int_id) {
  AKANTU_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <>
inline const InternalField<Real> &
Material::getInternal(const ID & int_id) const {
  auto it = internal_vectors_real.find(getID() + ":" + int_id);
  if (it == internal_vectors_real.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <>
inline InternalField<Real> & Material::getInternal(const ID & int_id) {
  auto it = internal_vectors_real.find(getID() + ":" + int_id);
  if (it == internal_vectors_real.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <>
inline const InternalField<UInt> &
Material::getInternal(const ID & int_id) const {
  auto it = internal_vectors_uint.find(getID() + ":" + int_id);
  if (it == internal_vectors_uint.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <>
inline InternalField<UInt> & Material::getInternal(const ID & int_id) {
  auto it = internal_vectors_uint.find(getID() + ":" + int_id);
  if (it == internal_vectors_uint.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const Array<T> & Material::getArray(const ID & vect_id, ElementType type,
                                           GhostType ghost_type) const {
  try {
    return this->template getInternal<T>(vect_id)(type, ghost_type);
  } catch (debug::Exception & e) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain a vector "
                                            << vect_id << " [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline Array<T> & Material::getArray(const ID & vect_id, ElementType type,
                                     GhostType ghost_type) {
  try {
    return this->template getInternal<T>(vect_id)(type, ghost_type);
  } catch (debug::Exception & e) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain a vector "
                                            << vect_id << " [" << e << "]");
  }
}

} // namespace akantu

#endif /* AKANTU_MATERIAL_INLINE_IMPL_HH_ */
