/**
 * @file   dumper_nodal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 26 2012
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Description of nodal fields
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
#ifndef AKANTU_DUMPER_NODAL_FIELD_HH_
#define AKANTU_DUMPER_NODAL_FIELD_HH_

/* -------------------------------------------------------------------------- */
#include "dumper_compute.hh"
#include "dumper_field.hh"
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumpers {

  /* ------------------------------------------------------------------------ */
  // This represents a iohelper compatible field
  template <typename T, bool filtered = false, class Container = Array<T>,
            class Filter = Array<UInt>>
  class NodalField : public dumpers::Field {
    /* ---------------------------------------------------------------------- */
    /* Typedefs                                                               */
    /* ---------------------------------------------------------------------- */
  public:
    using support_type = UInt;
    using types = TypeTraits<T, Vector<T>, Container>;

    class iterator : public iohelper::iterator<T, iterator, Vector<T>> {
    public:
      iterator(T * vect, UInt _offset, UInt _n, UInt _stride,
               const UInt * filter)
          : internal_it(vect), offset(_offset), n(_n), stride(_stride),
            filter(filter) {}

      bool operator!=(const iterator & it) const override {
        if (filter != nullptr) {
          return filter != it.filter;
        }
        return internal_it != it.internal_it;
      }

      iterator & operator++() override {
        if (filter != nullptr) {
          ++filter;
        } else {
          internal_it += offset;
        }
        return *this;
      }

      Vector<T> operator*() override {
        if (filter != nullptr) {
          return Vector<T>(internal_it + *(filter)*offset + stride, n);
        }
        return Vector<T>(internal_it + stride, n);
      }

    private:
      T * internal_it;
      UInt offset, n, stride;
      const UInt * filter{nullptr};
    };

    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  public:
    NodalField(const Container & _field, UInt _n = 0, UInt _stride = 0,
               const Filter * filter = nullptr)
        : field(_field), n(_n), stride(_stride), filter(filter), padding(0) {
      AKANTU_DEBUG_ASSERT(((not filtered) and filter == nullptr) or filtered,
                          "Filter passed to unfiltered NodalField!");
      AKANTU_DEBUG_ASSERT((filtered and this->filter != nullptr) or
                              (not filtered),
                          "No filter passed to filtered NodalField!");
      AKANTU_DEBUG_ASSERT(
          (filter != nullptr and this->filter->getNbComponent() == 1) or
              (filter == nullptr),
          "Multi-component filter given to NodalField ("
              << this->filter->getNbComponent()
              << " components detected, sould be 1");
      if (n == 0) {
        this->n = field.getNbComponent() - stride;
      }
    }

    /* ---------------------------------------------------------------------- */
    /* Methods                                                                */
    /* ---------------------------------------------------------------------- */
  public:
    void registerToDumper(const std::string & id,
                          iohelper::Dumper & dumper) override {
      dumper.addNodeDataField(id, *this);
    }

    inline iterator begin() {
      return iterator(field.storage(), field.getNbComponent(), n, stride,
                      filter == nullptr ? nullptr : filter->storage());
    }

    inline iterator end() {
      return iterator(field.storage() + field.getNbComponent() * field.size(),
                      field.getNbComponent(), n, stride,
                      filter == nullptr ? nullptr
                                        : filter->storage() + filter->size());
    }

    bool isHomogeneous() override { return true; }
    void checkHomogeneity() override { this->homogeneous = true; }

    virtual UInt getDim() {
      if (this->padding) {
        return this->padding;
      }
      return n;
    }

    void setPadding(UInt padding) { this->padding = padding; }

    UInt size() {
      if (filter != nullptr) {
        return filter->size();
      }
      return field.size();
    }

    inline std::shared_ptr<Field> connect(FieldComputeProxy & proxy) override {
      return proxy.connectToField(this);
    }

    /// for connection to a Homogenizer
    inline std::unique_ptr<ComputeFunctorInterface>
    connect(HomogenizerProxy & /*proxy*/) override {
      throw;
    }

    template <class T1 = T,
              std::enable_if_t<std::is_enum<T1>::value> * = nullptr>
    iohelper::DataType getDataType() {
      return iohelper::getDataType<UInt>();
    }

    template <class T1 = T,
              std::enable_if_t<not std::is_enum<T1>::value> * = nullptr>
    iohelper::DataType getDataType() {
      return iohelper::getDataType<T>();
    }

    /* ---------------------------------------------------------------------- */
    /* Class Members */
    /* ---------------------------------------------------------------------- */
  private:
    const Container & field;
    UInt n, stride;
    const Filter * filter{nullptr};
    UInt padding;
  };

} // namespace dumpers
} // namespace akantu
/* -------------------------------------------------------------------------- */
#endif /* AKANTU_DUMPER_NODAL_FIELD_HH_ */
