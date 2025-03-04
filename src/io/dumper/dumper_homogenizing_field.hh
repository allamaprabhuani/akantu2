/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_DUMPER_HOMOGENIZING_FIELD_HH_
#define AKANTU_DUMPER_HOMOGENIZING_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_compute.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumpers {

/* ------------------------------------------------------------------------ */
template <typename type>
inline type typeConverter(const type & /*input*/,
                          Vector<typename type::value_type> & /*res*/,
                          Int /*nb_data*/) {
  throw;
}

/* ------------------------------------------------------------------------ */
template <typename type>
inline Matrix<type> typeConverter(const Matrix<type> &input, Vector<type> &res,
                                  Int nb_data) {
  MatrixProxy<type> tmp(res.data(), input.rows(), nb_data / input.rows());
  return tmp;
}

/* ------------------------------------------------------------------------ */
template <typename type>
inline Vector<type> typeConverter(const Vector<type> & /*unused*/,
                                  Vector<type> &res, Int /*unused*/) {
  return res;
}

/* ------------------------------------------------------------------------ */
template <typename type>
class AvgHomogenizingFunctor : public ComputeFunctor<type, type> {
  /* ---------------------------------------------------------------------- */
  /* Typedefs */
  /* ---------------------------------------------------------------------- */
private:
  using value_type = typename type::value_type;

  /* ---------------------------------------------------------------------- */
  /* Constructors/Destructors */
  /* ---------------------------------------------------------------------- */
public:
  AvgHomogenizingFunctor(ElementTypeMap<Int> &nb_datas) {
    auto types = nb_datas.elementTypes();
    auto tit = types.begin();
    auto end = types.end();

    nb_data = nb_datas(*tit);

    for (; tit != end; ++tit) {
      if (nb_data != nb_datas(*tit)) {
        throw;
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /* Methods                                                                */
  /* ---------------------------------------------------------------------- */
public:
  type func(const type &d, Element /*global_index*/) override {
    Vector<value_type> res(this->nb_data);
    res.zero();
    if (d.size() % this->nb_data) {
      throw;
    }
    auto nb_to_average = d.size() / this->nb_data;

    auto &&ptr = d.data();
    for (Int i = 0; i < nb_to_average; ++i) {
      VectorProxy<const value_type> tmp(ptr, this->nb_data);
      res += tmp;
      ptr += this->nb_data;
    }
    res /= nb_to_average;
    return typeConverter(d, res, this->nb_data);
  };

  Int getDim() override { return nb_data; };
  Int getNbComponent(Int /*old_nb_comp*/) override { throw; };

  /* ---------------------------------------------------------------------- */
  /* Class Members */
  /* ---------------------------------------------------------------------- */

  /// The size of data: i.e. the size of the vector to be returned
  Int nb_data;
};
/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */
class HomogenizerProxy {
  /* ---------------------------------------------------------------------- */
  /* Constructors/Destructors */
  /* ---------------------------------------------------------------------- */
public:
  HomogenizerProxy() = default;

public:
  inline static std::unique_ptr<ComputeFunctorInterface>
  createHomogenizer(Field &field);

  template <typename T>
  inline std::unique_ptr<ComputeFunctorInterface> connectToField(T *field) {
    ElementTypeMap<Int> nb_components = field->getNbComponents();

    using ret_type = typename T::types::return_type;
    return this->instantiateHomogenizer<ret_type>(nb_components);
  }

  template <typename ret_type>
  inline std::unique_ptr<ComputeFunctorInterface>
  instantiateHomogenizer(ElementTypeMap<Int> &nb_components);
};

/* ------------------------------------------------------------------------ */

template <typename ret_type>
inline std::unique_ptr<ComputeFunctorInterface>
HomogenizerProxy::instantiateHomogenizer(ElementTypeMap<Int> &nb_components) {
  using Homogenizer = dumpers::AvgHomogenizingFunctor<ret_type>;
  return std::make_unique<Homogenizer>(nb_components);
}

template <>
inline std::unique_ptr<ComputeFunctorInterface>
HomogenizerProxy::instantiateHomogenizer<Vector<iohelper::ElemType>>(
    ElementTypeMap<Int> & /*nb_components*/) {
  throw;
  return nullptr;
}

/* ------------------------------------------------------------------------ */
/// for connection to a FieldCompute
template <typename SubFieldCompute, typename return_type,
          typename support_type_>
inline std::unique_ptr<ComputeFunctorInterface>
FieldCompute<SubFieldCompute, return_type, support_type_>::connect(
    HomogenizerProxy &proxy) {
  return proxy.connectToField(this);
}

template <typename SubFieldCompute, typename return_type>
inline std::unique_ptr<ComputeFunctorInterface>
FieldCompute<SubFieldCompute, return_type, Element>::connect(
    HomogenizerProxy &proxy) {
  return proxy.connectToField(this);
}

/* ------------------------------------------------------------------------ */
inline std::unique_ptr<ComputeFunctorInterface>
HomogenizerProxy::createHomogenizer(Field &field) {
  HomogenizerProxy homogenizer_proxy;
  return field.connect(homogenizer_proxy);
}

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_HOMOGENIZING_FIELD_HH_ */
