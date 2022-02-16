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
#include "integration_point.hh"
#include "material.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_MATERIAL_INLINE_IMPL_CC__
// #define __AKANTU_MATERIAL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline auto Material::addElement(ElementType type, Int element,
                                 GhostType ghost_type) {
  auto &el_filter = this->element_filter(type, ghost_type);
  el_filter.push_back(element);
  return el_filter.size() - 1;
}

/* -------------------------------------------------------------------------- */
inline auto Material::addElement(const Element &element) {
  return this->addElement(element.type, element.element, element.ghost_type);
}

/* -------------------------------------------------------------------------- */
constexpr inline Int Material::getCauchyStressMatrixSize(Int dim) {
  return (dim * dim);
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline void Material::gradUToF(const Eigen::MatrixBase<D1> &grad_u,
                                         Eigen::MatrixBase<D2> &F) {
  assert(F.size() >= grad_u.size() && grad_u.size() == dim * dim &&
         "The dimension of the tensor F should be greater or "
         "equal to the dimension of the tensor grad_u.");

  F.Identity();
  F.template block<dim, dim>(0, 0) += grad_u;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1>
constexpr inline decltype(auto)
Material::gradUToF(const Eigen::MatrixBase<D1> &grad_u) {
  Matrix<Real, dim, dim> F;
  gradUToF<dim>(grad_u, F);
  return F;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2, typename D3>
constexpr inline void Material::StoCauchy(const Eigen::MatrixBase<D1> &F,
                                          const Eigen::MatrixBase<D2> &S,
                                          Eigen::MatrixBase<D3> &sigma,
                                          const Real &C33) {
  Real J = F.determinant() * sqrt(C33);

  Matrix<Real, dim, dim> F_S;
  F_S = F * S;
  Real constant = J ? 1. / J : 0;
  sigma = constant * F_S * F.transpose();
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline decltype(auto)
Material::StoCauchy(const Eigen::MatrixBase<D1> &F,
                    const Eigen::MatrixBase<D2> &S, const Real &C33) {
  Matrix<Real, dim, dim> sigma;
  Material::StoCauchy<dim>(F, S, sigma, C33);
  return sigma;
}
/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
constexpr inline void Material::rightCauchy(const Eigen::MatrixBase<D1> &F,
                                            Eigen::MatrixBase<D2> &C) {
  C = F.transpose() * F;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D>
constexpr inline decltype(auto)
Material::rightCauchy(const Eigen::MatrixBase<D> &F) {
  Matrix<Real, dim, dim> C;
  rightCauchy(F, C);
  return C;
}

/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
constexpr inline void Material::leftCauchy(const Eigen::MatrixBase<D1> &F,
                                           Eigen::MatrixBase<D2> &B) {
  B = F * F.transpose();
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D>
constexpr inline decltype(auto)
Material::leftCauchy(const Eigen::MatrixBase<D> &F) {
  Matrix<Real, dim, dim> B;
  rightCauchy(F, B);
  return B;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline void
Material::gradUToEpsilon(const Eigen::MatrixBase<D1> &grad_u,
                         Eigen::MatrixBase<D2> &epsilon) {
  epsilon = .5 * (grad_u.transpose() + grad_u);
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1>
inline decltype(auto) constexpr Material::gradUToEpsilon(
    const Eigen::MatrixBase<D1> &grad_u) {
  Matrix<Real, dim, dim> epsilon;
  Material::gradUToEpsilon<dim>(grad_u, epsilon);
  return epsilon;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline void Material::gradUToE(const Eigen::MatrixBase<D1> &grad_u,
                                         Eigen::MatrixBase<D2> &E) {
  E = (grad_u.transpose() * grad_u + grad_u.transpose() + grad_u) / 2.;
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1>
constexpr inline decltype(auto)
Material::gradUToE(const Eigen::MatrixBase<D1> &grad_u) {
  Matrix<Real, dim, dim> E;
  gradUToE<dim>(grad_u, E);
  return E;
}

/* -------------------------------------------------------------------------- */
template <typename D1>
inline Real Material::stressToVonMises(const Eigen::MatrixBase<D1> &stress) {
  // compute deviatoric stress
  auto dim = stress.cols();
  auto &&deviatoric_stress =
      stress - Matrix<Real>::Identity(dim, dim) * stress.trace() / 3.;

  // return Von Mises stress
  return std::sqrt(3. * deviatoric_stress.doubleDot(deviatoric_stress) / 2.);
}

/* -------------------------------------------------------------------------- */
template <Int dim, typename D1, typename D2>
constexpr inline void
Material::setCauchyStressMatrix(const Eigen::MatrixBase<D1> &S_t,
                                Eigen::MatrixBase<D2> &sigma) {
  sigma.zero();

  /// see Finite ekement formulations for large deformation dynamic analysis,
  /// Bathe et al. IJNME vol 9, 1975, page 364 ^t \f$\tau\f$
  for (Int i = 0; i < dim; ++i) {
    for (Int m = 0; m < dim; ++m) {
      for (Int n = 0; n < dim; ++n) {
        sigma(i * dim + m, i * dim + n) = S_t(m, n);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline Element
Material::convertToLocalElement(const Element &global_element) const {
  auto ge = global_element.element;
#ifndef AKANTU_NDEBUG
  auto model_mat_index = this->model.getMaterialByElement(
      global_element.type, global_element.ghost_type)(ge);

  auto mat_index = this->model.getMaterialIndex(this->name);
  AKANTU_DEBUG_ASSERT(model_mat_index == mat_index,
                      "Conversion of a global  element in a local element for "
                      "the wrong material "
                          << this->name << std::endl);
#endif
  auto le = this->model.getMaterialLocalNumbering(
      global_element.type, global_element.ghost_type)(ge);

  Element tmp_quad{global_element.type, le, global_element.ghost_type};
  return tmp_quad;
}

/* -------------------------------------------------------------------------- */
inline Element
Material::convertToGlobalElement(const Element &local_element) const {
  auto le = local_element.element;
  auto ge =
      this->element_filter(local_element.type, local_element.ghost_type)(le);

  Element tmp_quad{local_element.type, ge, local_element.ghost_type};
  return tmp_quad;
}

/* -------------------------------------------------------------------------- */
inline IntegrationPoint
Material::convertToLocalPoint(const IntegrationPoint &global_point) const {
  const FEEngine &fem = this->model.getFEEngine();
  auto &&nb_quad = fem.getNbIntegrationPoints(global_point.type);
  auto &&el =
      this->convertToLocalElement(static_cast<const Element &>(global_point));
  return IntegrationPoint(el, global_point.num_point, nb_quad);
}

/* -------------------------------------------------------------------------- */
inline IntegrationPoint
Material::convertToGlobalPoint(const IntegrationPoint &local_point) const {
  const FEEngine &fem = this->model.getFEEngine();
  auto nb_quad = fem.getNbIntegrationPoints(local_point.type);
  Element el =
      this->convertToGlobalElement(static_cast<const Element &>(local_point));
  IntegrationPoint tmp_quad(el, local_point.num_point, nb_quad);
  return tmp_quad;
}

/* -------------------------------------------------------------------------- */
inline Int Material::getNbData(const Array<Element> &elements,
                               const SynchronizationTag &tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    return (this->isFiniteDeformation() ? 3 : 1) * spatial_dimension *
           spatial_dimension * sizeof(Real) *
           this->getModel().getNbIntegrationPoints(elements);
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void Material::packData(CommunicationBuffer &buffer,
                               const Array<Element> &elements,
                               const SynchronizationTag &tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      packElementDataHelper(piola_kirchhoff_2, buffer, elements);
      packElementDataHelper(gradu, buffer, elements);
    }
    packElementDataHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::unpackData(CommunicationBuffer &buffer,
                                 const Array<Element> &elements,
                                 const SynchronizationTag &tag) {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      unpackElementDataHelper(piola_kirchhoff_2, buffer, elements);
      unpackElementDataHelper(gradu, buffer, elements);
    }
    unpackElementDataHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline const Parameter &Material::getParam(const ID &param) const {
  try {
    return get(param);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the material "
                                     << getID());
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void Material::setParam(const ID &param, T value) {
  try {
    set<T>(param, value);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the material "
                                     << getID());
  }
  updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Material::packElementDataHelper(
    const ElementTypeMapArray<T> &data_to_pack, CommunicationBuffer &buffer,
    const Array<Element> &elements, const ID &fem_id) const {
  DataAccessor::packElementalDataHelper<T>(data_to_pack, buffer, elements, true,
                                           model.getFEEngine(fem_id));
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Material::unpackElementDataHelper(
    ElementTypeMapArray<T> &data_to_unpack, CommunicationBuffer &buffer,
    const Array<Element> &elements, const ID &fem_id) {
  DataAccessor::unpackElementalDataHelper<T>(data_to_unpack, buffer, elements,
                                             true, model.getFEEngine(fem_id));
}

/* -------------------------------------------------------------------------- */
template <>
inline void Material::registerInternal<Real>(InternalField<Real> &vect) {
  internal_vectors_real[vect.getID()] = &vect;
}

template <>
inline void Material::registerInternal<Int>(InternalField<Int> &vect) {
  internal_vectors_int[vect.getID()] = &vect;
}

template <>
inline void Material::registerInternal<bool>(InternalField<bool> &vect) {
  internal_vectors_bool[vect.getID()] = &vect;
}

/* -------------------------------------------------------------------------- */
template <>
inline void Material::unregisterInternal<Real>(InternalField<Real> &vect) {
  internal_vectors_real.erase(vect.getID());
}

template <>
inline void Material::unregisterInternal<Int>(InternalField<Int> &vect) {
  internal_vectors_int.erase(vect.getID());
}

template <>
inline void Material::unregisterInternal<bool>(InternalField<bool> &vect) {
  internal_vectors_bool.erase(vect.getID());
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline bool Material::isInternal(const ID & /*id*/,
                                 ElementKind /*element_kind*/) const {
  AKANTU_TO_IMPLEMENT();
}

template <>
inline bool Material::isInternal<Real>(const ID &id,
                                       ElementKind element_kind) const {
  auto internal_array = internal_vectors_real.find(this->getID() + ":" + id);

  return not(internal_array == internal_vectors_real.end() ||
             internal_array->second->getElementKind() != element_kind);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline ElementTypeMap<Int>
Material::getInternalDataPerElem(const ID &field_id,
                                 ElementKind element_kind) const {

  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in material "
                                                   << this->name);
  }

  const InternalField<T> &internal_field =
      this->template getInternal<T>(field_id);
  const auto &fe_engine = internal_field.getFEEngine();
  auto nb_data_per_quad = internal_field.getNbComponent();

  ElementTypeMap<Int> res;
  for (auto ghost_type : ghost_types) {
    for (auto &&type : internal_field.elementTypes(ghost_type)) {
      auto nb_quadrature_points =
          fe_engine.getNbIntegrationPoints(type, ghost_type);
      res(type, ghost_type) = nb_data_per_quad * nb_quadrature_points;
    }
  }

  return res;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Material::flattenInternal(const std::string &field_id,
                               ElementTypeMapArray<T> &internal_flat,
                               const GhostType ghost_type,
                               ElementKind element_kind) const {

  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in material "
                                                   << this->name);
  }

  const auto &internal_field = this->template getInternal<T>(field_id);

  const auto &fe_engine = internal_field.getFEEngine();
  const auto &mesh = fe_engine.getMesh();

  for (auto &&type : internal_field.filterTypes(ghost_type)) {
    const auto &src_vect = internal_field(type, ghost_type);
    const auto &filter = internal_field.getFilter(type, ghost_type);

    // total number of elements in the corresponding mesh
    auto nb_element_dst = mesh.getNbElement(type, ghost_type);
    // number of element in the internal field
    auto nb_element_src = filter.size();
    // number of quadrature points per elem
    auto nb_quad_per_elem = fe_engine.getNbIntegrationPoints(type);
    // number of data per quadrature point
    auto nb_data_per_quad = internal_field.getNbComponent();

    if (not internal_flat.exists(type, ghost_type)) {
      internal_flat.alloc(nb_element_dst * nb_quad_per_elem, nb_data_per_quad,
                          type, ghost_type);
    }

    if (nb_element_src == 0) {
      continue;
    }

    // number of data per element
    auto nb_data = nb_quad_per_elem * nb_data_per_quad;

    auto &dst_vect = internal_flat(type, ghost_type);
    dst_vect.resize(nb_element_dst * nb_quad_per_elem);

    auto it_dst = make_view(dst_vect, nb_data).begin();

    for (auto &&data : zip(filter, make_view(src_vect, nb_data))) {
      it_dst[std::get<0>(data)] = std::get<1>(data);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const InternalField<T> &
Material::getInternal(const ID & /*int_id*/) const {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline InternalField<T> &Material::getInternal(const ID & /*int_id*/) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline const InternalField<Real> &
Material::getInternal(const ID &int_id) const {
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
inline InternalField<Real> &Material::getInternal(const ID &int_id) {
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
inline const InternalField<Int> &Material::getInternal(const ID &int_id) const {
  auto it = internal_vectors_int.find(getID() + ":" + int_id);
  if (it == internal_vectors_int.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <> inline InternalField<Int> &Material::getInternal(const ID &int_id) {
  auto it = internal_vectors_int.find(getID() + ":" + int_id);
  if (it == internal_vectors_int.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const Array<T> &Material::getArray(const ID &vect_id, ElementType type,
                                          GhostType ghost_type) const {
  try {
    return this->template getInternal<T>(vect_id)(type, ghost_type);
  } catch (debug::Exception &e) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain a vector "
                                            << vect_id << " [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline Array<T> &Material::getArray(const ID &vect_id, ElementType type,
                                    GhostType ghost_type) {
  try {
    return this->template getInternal<T>(vect_id)(type, ghost_type);
  } catch (debug::Exception &e) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain a vector "
                                            << vect_id << " [" << e << "]");
  }
}

} // namespace akantu

//#endif /* __AKANTU_MATERIAL_INLINE_IMPL_CC__ */
