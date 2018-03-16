/**
 * @file   data_accessor.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Sun Feb 04 2018
 *
 * @brief  Interface of accessors for pack_unpack system
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
#include "aka_common.hh"
#include "communication_buffer.hh"
#include "element.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DATA_ACCESSOR_HH__
#define __AKANTU_DATA_ACCESSOR_HH__

namespace akantu {
class FEEngine;
} // akantu

namespace akantu {

class DataAccessorBase {
public:
  DataAccessorBase() = default;
  virtual ~DataAccessorBase() = default;
};

template <class T> class DataAccessor : public virtual DataAccessorBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DataAccessor() = default;
  ~DataAccessor() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * @brief get  the number of  data to exchange  for a given array of T
   * (elements or dofs) and a given akantu::SynchronizationTag
   */
  virtual UInt getNbData(const Array<T> & elements,
                         const SynchronizationTag & tag) const = 0;

  /**
   * @brief pack the data for a given array of T (elements or dofs) and a given
   * akantu::SynchronizationTag
   */
  virtual void packData(CommunicationBuffer & buffer, const Array<T> & element,
                        const SynchronizationTag & tag) const = 0;

  /**
   * @brief unpack the data for a given array of T (elements or dofs) and a
   * given akantu::SynchronizationTag
   */
  virtual void unpackData(CommunicationBuffer & buffer,
                          const Array<T> & element,
                          const SynchronizationTag & tag) = 0;
};

/* -------------------------------------------------------------------------- */
/* Specialization                                                             */
/* -------------------------------------------------------------------------- */
template <> class DataAccessor<Element> : public virtual DataAccessorBase {
public:
  DataAccessor() = default;
  ~DataAccessor() override = default;

  virtual UInt getNbData(const Array<Element> & elements,
                         const SynchronizationTag & tag) const = 0;
  virtual void packData(CommunicationBuffer & buffer,
                        const Array<Element> & element,
                        const SynchronizationTag & tag) const = 0;
  virtual void unpackData(CommunicationBuffer & buffer,
                          const Array<Element> & element,
                          const SynchronizationTag & tag) = 0;

  /* ------------------------------------------------------------------------ */
public:
  template <typename T, bool pack_helper>
  static void
  packUnpackNodalDataHelper(Array<T> & data, CommunicationBuffer & buffer,
                            const Array<Element> & elements, const Mesh & mesh);

  /* ------------------------------------------------------------------------ */
  template <typename T, bool pack_helper>
  static void packUnpackElementalDataHelper(
      ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
      const Array<Element> & element, bool per_quadrature_point_data,
      const FEEngine & fem);

  /* ------------------------------------------------------------------------ */
  template <typename T>
  static void
  packNodalDataHelper(const Array<T> & data, CommunicationBuffer & buffer,
                      const Array<Element> & elements, const Mesh & mesh) {
    packUnpackNodalDataHelper<T, true>(const_cast<Array<T> &>(data), buffer,
                                       elements, mesh);
  }

  template <typename T>
  static inline void
  unpackNodalDataHelper(Array<T> & data, CommunicationBuffer & buffer,
                        const Array<Element> & elements, const Mesh & mesh) {
    packUnpackNodalDataHelper<T, false>(data, buffer, elements, mesh);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  static inline void
  packElementalDataHelper(const ElementTypeMapArray<T> & data_to_pack,
                          CommunicationBuffer & buffer,
                          const Array<Element> & elements,
                          bool per_quadrature_point, const FEEngine & fem) {
    packUnpackElementalDataHelper<T, true>(
        const_cast<ElementTypeMapArray<T> &>(data_to_pack), buffer, elements,
        per_quadrature_point, fem);
  }

  template <typename T>
  static inline void
  unpackElementalDataHelper(ElementTypeMapArray<T> & data_to_unpack,
                            CommunicationBuffer & buffer,
                            const Array<Element> & elements,
                            bool per_quadrature_point, const FEEngine & fem) {
    packUnpackElementalDataHelper<T, false>(data_to_unpack, buffer, elements,
                                            per_quadrature_point, fem);
  }
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <> class DataAccessor<UInt> : public virtual DataAccessorBase {
public:
  DataAccessor() = default;
  ~DataAccessor() override = default;

  virtual UInt getNbData(const Array<UInt> & elements,
                         const SynchronizationTag & tag) const = 0;
  virtual void packData(CommunicationBuffer & buffer,
                        const Array<UInt> & element,
                        const SynchronizationTag & tag) const = 0;
  virtual void unpackData(CommunicationBuffer & buffer,
                          const Array<UInt> & element,
                          const SynchronizationTag & tag) = 0;
  /* ------------------------------------------------------------------------ */
public:
  template <typename T, bool pack_helper>
  static void packUnpackDOFDataHelper(Array<T> & data,
                                      CommunicationBuffer & buffer,
                                      const Array<UInt> & dofs);

  template <typename T>
  static inline void packDOFDataHelper(const Array<T> & data_to_pack,
                                       CommunicationBuffer & buffer,
                                       const Array<UInt> & dofs) {
    packUnpackDOFDataHelper<T, true>(const_cast<Array<T> &>(data_to_pack),
                                     buffer, dofs);
  }

  template <typename T>
  static inline void unpackDOFDataHelper(Array<T> & data_to_unpack,
                                         CommunicationBuffer & buffer,
                                         const Array<UInt> & dofs) {
    packUnpackDOFDataHelper<T, false>(data_to_unpack, buffer, dofs);
  }
};

/* -------------------------------------------------------------------------- */
template <typename T> class AddOperation {
public:
  inline T operator()(T & a, T & b) { return a + b; };
};

template <typename T> class IdentityOperation {
public:
  inline T & operator()(T &, T & b) { return b; };
};
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <template <class> class Op, class T>
class ReduceUIntDataAccessor : public virtual DataAccessor<UInt> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ReduceUIntDataAccessor(Array<T> & data, const SynchronizationTag & tag)
      : data(data), tag(tag) {}

  ~ReduceUIntDataAccessor() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  UInt getNbData(const Array<UInt> & elements,
                 const SynchronizationTag & tag) const override {
    if (tag != this->tag)
      return 0;

    Vector<T> tmp(data.getNbComponent());
    return elements.size() * CommunicationBuffer::sizeInBuffer(tmp);
  }

  /* ------------------------------------------------------------------------ */
  void packData(CommunicationBuffer & buffer, const Array<UInt> & elements,
                const SynchronizationTag & tag) const override {
    if (tag != this->tag)
      return;

    auto data_it = data.begin(data.getNbComponent());
    for (auto el : elements) {
      buffer << Vector<T>(data_it[el]);
    }
  }

  /* ------------------------------------------------------------------------ */
  void unpackData(CommunicationBuffer & buffer, const Array<UInt> & elements,
                  const SynchronizationTag & tag) override {
    if (tag != this->tag)
      return;

    auto data_it = data.begin(data.getNbComponent());
    for (auto el : elements) {
      Vector<T> unpacked(data.getNbComponent());
      Vector<T> vect(data_it[el]);
      buffer >> unpacked;

      vect = oper(vect, unpacked);
    }
  }

protected:
  /// data to (un)pack
  Array<T> & data;

  /// Tag to consider
  SynchronizationTag tag;

  /// reduction operator
  Op<Vector<T>> oper;
};

/* -------------------------------------------------------------------------- */
template <class T>
using SimpleUIntDataAccessor = ReduceUIntDataAccessor<IdentityOperation, T>;

} // akantu

#endif /* __AKANTU_DATA_ACCESSOR_HH__ */
