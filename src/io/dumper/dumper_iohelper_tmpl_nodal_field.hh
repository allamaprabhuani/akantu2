/**
 * @file   dumper_iohelper_tmpl_nodal_field.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Oct 25 14:16:35 2012
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
template< typename T, template<typename> class return_type>
class DumperIOHelper::NodalField : public Field {
public:
  /* -----------------------------------------------------------------------*/
  class iterator : iohelper::iterator< T, iterator, return_type<T> > {
  protected:
    typedef typename Vector<T>::template const_iterator< return_type<T> > internal_iterator;
  public:
    iterator(const internal_iterator & it) : vit(it) {}

    bool           operator!=(const iterator & it) const { return vit != it.vit; }
    iterator &     operator++() { ++vit; return *this; };
    return_type<T> operator* (){ return *vit; };
  private:
    internal_iterator vit;
  };

  /* ---------------------------------------------------------------------- */
  NodalField(const Vector<T> & field, UInt n = 0, UInt m = 1) : field(field), n(n), m(m) {
    if(n == 0) { this->n = field.getNbComponent(); }
  }

  virtual void registerToDumper(const std::string & id,
				iohelper::Dumper & dumper) {
    dumper.addNodeDataField(id, *this);
  }


  inline iterator begin() { return iterator_helper<T, return_type>::begin(field, n ,m, field.getSize()); }
  inline iterator end  () { return iterator_helper<T, return_type>::end(field, n ,m, field.getSize()); }

  bool isHomogeneous() { return true; }

  virtual UInt getDim() {
    if(padding_n && padding_m)
      return padding_m*padding_n;
    else return n*m;
  }

  UInt size() { return field.getSize(); }
private:
  const Vector<T> & field;
  UInt n, m;
};


template<>
inline typename DumperIOHelper::NodalField<Real, types::Matrix>::iterator
DumperIOHelper::NodalField<Real, types::Matrix>::begin() {
  return iterator(field.begin(n, m));
}

template<>
inline typename DumperIOHelper::NodalField<Real, types::Matrix>::iterator
DumperIOHelper::NodalField<Real, types::Matrix>::end() {
  return iterator(field.end(n, m));
}
