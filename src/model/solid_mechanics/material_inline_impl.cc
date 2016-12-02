/**
 * @file   material_inline_impl.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Wed Nov 25 2015
 *
 * @brief  Implementation of the inline functions of the class material
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_INLINE_IMPL_CC__
#define __AKANTU_MATERIAL_INLINE_IMPL_CC__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline UInt Material::addElement(const ElementType & type, UInt element,
                                 const GhostType & ghost_type) {
  Array<UInt> & el_filter = this->element_filter(type, ghost_type);
  el_filter.push_back(element);
  return el_filter.getSize() - 1;
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getTangentStiffnessVoigtSize(UInt dim) const {
  return (dim * (dim - 1) / 2 + dim);
}
/* -------------------------------------------------------------------------- */
inline UInt Material::getCauchyStressMatrixSize(UInt dim) const {
  return (dim * dim);
}
/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::gradUToF(const Matrix<Real> & grad_u, Matrix<Real> & F) {
  AKANTU_DEBUG_ASSERT(F.size() >= grad_u.size() && grad_u.size() == dim * dim,
                      "The dimension of the tensor F should be greater or "
                      "equal to the dimension of the tensor grad_u.");

  F.eye();

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      F(i, j) += grad_u(i, j);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::computeCauchyStressOnQuad(const Matrix<Real> & F,
                                                const Matrix<Real> & piola,
                                                Matrix<Real> & sigma,
                                                const Real & C33) const {

  Real J = F.det() * sqrt(C33);

  Matrix<Real> F_S(dim, dim);
  F_S.mul<false, false>(F, piola);
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
  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      epsilon(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::gradUToGreenStrain(const Matrix<Real> & grad_u,
                                         Matrix<Real> & epsilon) {
  epsilon.mul<true, false>(grad_u, grad_u, .5);

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      epsilon(i, j) += 0.5 * (grad_u(i, j) + grad_u(j, i));
}

/* -------------------------------------------------------------------------- */
inline Real Material::stressToVonMises(const Matrix<Real> & stress) {
  // compute deviatoric stress
  UInt dim = stress.cols();
  Matrix<Real> deviatoric_stress =
      Matrix<Real>::eye(dim, -1. * stress.trace() / 3.);

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      deviatoric_stress(i, j) += stress(i, j);

  // return Von Mises stress
  return std::sqrt(3. * deviatoric_stress.doubleDot(deviatoric_stress) / 2.);
}

/* ---------------------------------------------------------------------------*/
template <UInt dim>
inline void Material::setCauchyStressArray(const Matrix<Real> & S_t,
                                           Matrix<Real> & sigma_voight) {
  AKANTU_DEBUG_IN();
  sigma_voight.clear();
  // see Finite ekement formulations for large deformation dynamic analysis,
  // Bathe et al. IJNME vol 9, 1975, page 364 ^t\tau

  /*
   * 1d: [ s11 ]'
   * 2d: [ s11 s22 s12 ]'
   * 3d: [ s11 s22 s33 s23 s13 s12 ]
   */
  for (UInt i = 0; i < dim; ++i) // diagonal terms
    sigma_voight(i, 0) = S_t(i, i);

  for (UInt i = 1; i < dim; ++i) // term s12 in 2D and terms s23 s13 in 3D
    sigma_voight(dim + i - 1, 0) = S_t(dim - i - 1, dim - 1);

  for (UInt i = 2; i < dim; ++i) // term s13 in 3D
    sigma_voight(dim + i, 0) = S_t(0, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void Material::setCauchyStressMatrix(const Matrix<Real> & S_t,
                                            Matrix<Real> & sigma) {
  AKANTU_DEBUG_IN();

  sigma.clear();

  /// see Finite ekement formulations for large deformation dynamic analysis,
  /// Bathe et al. IJNME vol 9, 1975, page 364 ^t\tau
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
inline Element
Material::convertToLocalElement(const Element & global_element) const {
  UInt ge = global_element.element;
#ifndef AKANTU_NDEBUG
  UInt model_mat_index = this->model->getMaterialByElement(
      global_element.type, global_element.ghost_type)(ge);

  UInt mat_index = this->model->getMaterialIndex(this->name);
  AKANTU_DEBUG_ASSERT(model_mat_index == mat_index,
                      "Conversion of a global  element in a local element for "
                      "the wrong material "
                          << this->name << std::endl);
#endif
  UInt le = this->model->getMaterialLocalNumbering(
      global_element.type, global_element.ghost_type)(ge);

  Element tmp_quad(global_element.type, le, global_element.ghost_type,
                   global_element.kind);
  return tmp_quad;
}

/* -------------------------------------------------------------------------- */
inline Element
Material::convertToGlobalElement(const Element & local_element) const {
  UInt le = local_element.element;
  UInt ge =
      this->element_filter(local_element.type, local_element.ghost_type)(le);

  Element tmp_quad(local_element.type, ge, local_element.ghost_type,
                   local_element.kind);
  return tmp_quad;
}

/* -------------------------------------------------------------------------- */
inline IntegrationPoint
Material::convertToLocalPoint(const IntegrationPoint & global_point) const {
  const FEEngine & fem = this->model->getFEEngine();
  UInt nb_quad = fem.getNbIntegrationPoints(global_point.type);
  Element el =
      this->convertToLocalElement(static_cast<const Element &>(global_point));
  IntegrationPoint tmp_quad(el, global_point.num_point, nb_quad);
  return tmp_quad;
}

/* -------------------------------------------------------------------------- */
inline IntegrationPoint
Material::convertToGlobalPoint(const IntegrationPoint & local_point) const {
  const FEEngine & fem = this->model->getFEEngine();
  UInt nb_quad = fem.getNbIntegrationPoints(local_point.type);
  Element el =
      this->convertToGlobalElement(static_cast<const Element &>(local_point));
  IntegrationPoint tmp_quad(el, local_point.num_point, nb_quad);
  return tmp_quad;
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getNbData(const Array<Element> & elements,
                                const SynchronizationTag & tag) const {
  if (tag == _gst_smm_stress) {
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
  if (tag == _gst_smm_stress) {
    if (this->isFiniteDeformation()) {
      packElementDataHelper(piola_kirchhoff_2, buffer, elements);
      packElementDataHelper(gradu, buffer, elements);
    }
    packElementDataHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::unpackData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const SynchronizationTag & tag) {
  if (tag == _gst_smm_stress) {
    if (this->isFiniteDeformation()) {
      unpackElementDataHelper(piola_kirchhoff_2, buffer, elements);
      unpackElementDataHelper(gradu, buffer, elements);
    }
    unpackElementDataHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline const Parameter & Material::getParam(const ID & param) const {
  try {
    return get(param);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the material "
                                     << getID());
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Material::setParam(const ID & param, T value) {
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
    const ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & elements, const ID & fem_id) const {
  DataAccessor::packElementalDataHelper<T>(data_to_pack, buffer, elements, true,
                                           model->getFEEngine(fem_id));
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Material::unpackElementDataHelper(
    ElementTypeMapArray<T> & data_to_unpack, CommunicationBuffer & buffer,
    const Array<Element> & elements, const ID & fem_id) {
  DataAccessor::unpackElementalDataHelper<T>(data_to_unpack, buffer, elements,
                                             true, model->getFEEngine(fem_id));
}

/* -------------------------------------------------------------------------- */
template <>
inline void Material::registerInternal<Real>(InternalField<Real> & vect) {
  internal_vectors_real[vect.getID()] = &vect;
}

template <>
inline void Material::registerInternal<UInt>(InternalField<UInt> & vect) {
  internal_vectors_uint[vect.getID()] = &vect;
}

template <>
inline void Material::registerInternal<bool>(InternalField<bool> & vect) {
  internal_vectors_bool[vect.getID()] = &vect;
}

/* -------------------------------------------------------------------------- */
template <>
inline void Material::unregisterInternal<Real>(InternalField<Real> & vect) {
  internal_vectors_real.erase(vect.getID());
}

template <>
inline void Material::unregisterInternal<UInt>(InternalField<UInt> & vect) {
  internal_vectors_uint.erase(vect.getID());
}

template <>
inline void Material::unregisterInternal<bool>(InternalField<bool> & vect) {
  internal_vectors_bool.erase(vect.getID());
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline bool Material::isInternal(__attribute__((unused)) const ID & id,
                                 __attribute__((unused)) const ElementKind & element_kind) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

template <>
inline bool Material::isInternal<Real>(const ID & id,
                                       const ElementKind & element_kind) const {
  std::map<ID, InternalField<Real> *>::const_iterator internal_array =
      internal_vectors_real.find(this->getID() + ":" + id);

  if (internal_array == internal_vectors_real.end() ||
      internal_array->second->getElementKind() != element_kind)
    return false;
  return true;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline ElementTypeMap<UInt>
Material::getInternalDataPerElem(const ID & field_id,
                                 const ElementKind & element_kind) const {

  if (!this->template isInternal<T>(field_id, element_kind))
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in material "
                                                   << this->name);

  const InternalField<T> & internal_field =
      this->template getInternal<T>(field_id);
  const FEEngine & fe_engine = internal_field.getFEEngine();
  UInt nb_data_per_quad = internal_field.getNbComponent();

  ElementTypeMap<UInt> res;
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    typedef typename InternalField<T>::type_iterator type_iterator;
    type_iterator tit = internal_field.firstType(*gt);
    type_iterator tend = internal_field.lastType(*gt);

    for (; tit != tend; ++tit) {
      UInt nb_quadrature_points = fe_engine.getNbIntegrationPoints(*tit, *gt);
      res(*tit, *gt) = nb_data_per_quad * nb_quadrature_points;
    }
  }
  return res;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Material::flattenInternal(const std::string & field_id,
                               ElementTypeMapArray<T> & internal_flat,
                               const GhostType ghost_type,
                               ElementKind element_kind) const {

  if (!this->template isInternal<T>(field_id, element_kind))
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in material "
                                                   << this->name);

  const InternalField<T> & internal_field =
      this->template getInternal<T>(field_id);

  const FEEngine & fe_engine = internal_field.getFEEngine();
  const Mesh & mesh = fe_engine.getMesh();

  typedef typename InternalField<T>::filter_type_iterator type_iterator;
  type_iterator tit = internal_field.filterFirstType(ghost_type);
  type_iterator tend = internal_field.filterLastType(ghost_type);

  for (; tit != tend; ++tit) {
    ElementType type = *tit;

    const Array<Real> & src_vect = internal_field(type, ghost_type);
    const Array<UInt> & filter = internal_field.getFilter(type, ghost_type);

    // total number of elements in the corresponding mesh
    UInt nb_element_dst = mesh.getNbElement(type, ghost_type);
    // number of element in the internal field
    UInt nb_element_src = filter.getSize();
    // number of quadrature points per elem
    UInt nb_quad_per_elem = fe_engine.getNbIntegrationPoints(type);
    // number of data per quadrature point
    UInt nb_data_per_quad = internal_field.getNbComponent();

    if (!internal_flat.exists(type, ghost_type)) {
      internal_flat.alloc(nb_element_dst * nb_quad_per_elem, nb_data_per_quad,
                          type, ghost_type);
    }

    if (nb_element_src == 0)
      continue;

    // number of data per element
    UInt nb_data = nb_quad_per_elem * nb_data_per_quad;

    Array<Real> & dst_vect = internal_flat(type, ghost_type);
    dst_vect.resize(nb_element_dst * nb_quad_per_elem);

    Array<UInt>::const_scalar_iterator it = filter.begin();
    Array<UInt>::const_scalar_iterator end = filter.end();

    Array<Real>::const_vector_iterator it_src =
        src_vect.begin_reinterpret(nb_data, nb_element_src);
    Array<Real>::vector_iterator it_dst =
        dst_vect.begin_reinterpret(nb_data, nb_element_dst);

    for (; it != end; ++it, ++it_src) {
      it_dst[*it] = *it_src;
    }
  }
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_INLINE_IMPL_CC__ */
