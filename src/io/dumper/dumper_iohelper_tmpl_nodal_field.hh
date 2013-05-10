/**
 * @file   dumper_iohelper_tmpl_nodal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 26 21:52:40 2012
 *
 * @brief  Description of nodal fields
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
template<typename T, class Container, class Filter>
class DumperIOHelper::NodalField<T, false, Container, Filter> : public Field {
public:
  /* -----------------------------------------------------------------------*/
  class iterator : public iohelper::iterator< T, iterator, Vector<T> > {
  public:
    iterator(T * vect, UInt offset, UInt n, UInt stride, __attribute__ ((unused)) const UInt * filter = NULL) :
      internal_it(vect), offset(offset), n(n), stride(stride) {}

    bool operator!=(const iterator & it) const { return internal_it != it.internal_it; }
    iterator & operator++() { internal_it += offset; return *this; };
    Vector<T> operator* (){ return Vector<T>(internal_it + stride, n); };
  private:
    T * internal_it;
    UInt offset, n, stride;
  };

  /* ---------------------------------------------------------------------- */
  NodalField(const Container & field, UInt n = 0, UInt stride = 0, __attribute__ ((unused)) const Filter * filter = NULL) :
    field(field), n(n), stride(stride) {
    AKANTU_DEBUG_ASSERT(filter == NULL, "Filter passed to unfiltered NodalField!");
    if(n == 0) { this->n = field.getNbComponent() - stride; }
  }

  virtual void registerToDumper(const std::string & id,
				iohelper::Dumper & dumper) {
    dumper.addNodeDataField(id, *this);
  }

  inline iterator begin() {
    return iterator(field.storage(), field.getNbComponent(), n, stride);
  }

  inline iterator end  () {
    return iterator(field.storage() + field.getNbComponent()*field.getSize(),
		    field.getNbComponent(), n, stride);
  }

  bool isHomogeneous() { return true; }

  virtual UInt getDim() {
    if(this->padding_n && this->padding_m)
      return this->padding_n * this->padding_m;
    else return n;
  }

  UInt size() { return field.getSize(); }

  iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }

private:
  const Container & field;
  UInt n, stride;
};



/* -------------------------------------------------------------------------- */
template<typename T, class Container, class Filter>
class DumperIOHelper::NodalField<T, true, Container, Filter> : public Field {
public:
  /* ------------------------------------------------------------------------ */
  class iterator : public iohelper::iterator< T, iterator, Vector<T> > {

  public:
    iterator(T * const vect, UInt _offset, UInt _n, UInt _stride, const UInt * filter) :
      internal_it(vect), offset(_offset), n(_n), stride(_stride),
      filter(filter) {}

    bool operator!=(const iterator & it) const {
      return filter != it.filter;
    }

    iterator & operator++() {
      ++filter;
      return *this;
    }

    Vector<T> operator* () {
      return Vector<T>(internal_it + *(filter)*offset + stride, n);
    }
  private:
    T * const internal_it;
    UInt offset, n, stride;
    const UInt * filter;
  };

  /* ---------------------------------------------------------------------- */
  NodalField(const Container & _field,
             UInt _n = 0, UInt _stride = 0,
             const Filter * filter = NULL)
  : field(_field), n(_n), stride(_stride), filter(filter) {
    AKANTU_DEBUG_ASSERT(filter != NULL, "No filter passed to filtered NodalField!");
    AKANTU_DEBUG_ASSERT(filter->getNbComponent()==1, "Multi-component filter given to NodalField!");
    if(n == 0) {
      this->n = field.getNbComponent() - stride;
    }
  }

  virtual void registerToDumper(const std::string & id, iohelper::Dumper & dumper) {
    dumper.addNodeDataField(id, *this);
  }

  inline iterator begin() {
    return iterator(field.storage(), field.getNbComponent(), n, stride, filter->storage());
  }

  inline iterator end() {
    return iterator(field.storage(), field.getNbComponent(), n, stride, filter->storage()+filter->getSize());
  }

  bool isHomogeneous() {
    return true;
  }

  virtual UInt getDim() {
    if(this->padding_n && this->padding_m)
      return this->padding_n * this->padding_m;
    else return n;
  }

  UInt size() {
    return filter->getSize();
    }

  iohelper::DataType getDataType() {
    return iohelper::getDataType<T>();
  }

private:
  const Container & field;
  UInt n, stride;
  const Filter * filter;

};
