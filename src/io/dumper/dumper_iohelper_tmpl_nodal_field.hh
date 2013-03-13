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
template<typename T>
class DumperIOHelper::NodalField : public Field {
public:
  /* -----------------------------------------------------------------------*/
  class iterator : public iohelper::iterator< T, iterator, Vector<T> > {
  protected:
    typedef typename Array<T>::template const_iterator< Vector<T> > internal_iterator;
  public:
    iterator(T * vect, UInt offset, UInt n, UInt stride) :
      internal_it(vect), offset(offset), n(n), stride(stride) {}

    bool operator!=(const iterator & it) const { return internal_it != it.internal_it; }
    iterator & operator++() { internal_it += offset; return *this; };
    Vector<T> operator* (){ return Vector<T>(internal_it + stride, n); };
  private:
    T * internal_it;
    UInt offset, n, stride;
  };

  /* ---------------------------------------------------------------------- */
  NodalField(const Array<T> & field, UInt n = 0, UInt stride = 0) :
    field(field), n(n), stride(stride) {
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
  const Array<T> & field;
  UInt n, stride;
};

