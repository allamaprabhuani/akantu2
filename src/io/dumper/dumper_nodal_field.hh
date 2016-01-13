/**
 * @file   dumper_nodal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  Description of nodal fields
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_DUMPER_NODAL_FIELD_HH__
#define __AKANTU_DUMPER_NODAL_FIELD_HH__

#include "dumper_field.hh"
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__

// This represents a iohelper compatible field
template<typename T, bool filtered = false,
	 class Container = Array<T>, class Filter = Array<UInt> >
class NodalField;

/* -------------------------------------------------------------------------- */
template<typename T, class Container, class Filter>
class NodalField<T, false, Container, Filter> : public dumper::Field {
public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
 
  /// associated iterator with any nodal field (non filetered)
  class iterator : public iohelper::iterator< T, iterator, Vector<T> > {
  public:
    iterator(T * vect, UInt offset, UInt n, UInt stride,
	     __attribute__ ((unused)) const UInt * filter = NULL) :

      internal_it(vect), offset(offset), n(n), stride(stride) {}

    bool operator!=(const iterator & it) const { return internal_it != it.internal_it; }
    iterator & operator++() { internal_it += offset; return *this; };
    Vector<T> operator* (){ return Vector<T>(internal_it + stride, n); };
  private:
    T * internal_it;
    UInt offset, n, stride;
  };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  NodalField(const Container & field, UInt n = 0, UInt stride = 0,
	     __attribute__ ((unused)) const Filter * filter = NULL) :

    field(field), n(n), stride(stride), padding(0) {
    AKANTU_DEBUG_ASSERT(filter == NULL, "Filter passed to unfiltered NodalField!");
    if(n == 0) { this->n = field.getNbComponent() - stride; }
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  
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
  void checkHomogeneity() { this->homogeneous = true; }

  virtual UInt getDim() {
    if(this->padding) return this->padding;
    else              return n;
  }

  void setPadding(UInt padding){this->padding = padding;}

  UInt size() { return field.getSize(); }

  iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  

private:
  const Container & field;
  UInt n, stride;
  UInt padding;
};



/* -------------------------------------------------------------------------- */
template<typename T, class Container, class Filter>
class NodalField<T, true, Container, Filter> : public dumper::Field {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
  
public:
  class iterator : public iohelper::iterator< T, iterator, Vector<T> > {

  public:
    iterator(T * const vect, UInt _offset, UInt _n, 
	     UInt _stride, const UInt * filter) :

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

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  NodalField(const Container & _field,
             UInt _n = 0, UInt _stride = 0,
             const Filter * filter = NULL)
    : field(_field), n(_n), stride(_stride), filter(filter), padding(0) {
    AKANTU_DEBUG_ASSERT(this->filter != NULL, 
			"No filter passed to filtered NodalField!");

    AKANTU_DEBUG_ASSERT(this->filter->getNbComponent()==1, 
			"Multi-component filter given to NodalField (" 
			<< this->filter->getNbComponent() 
			<< " components detected, sould be 1");
    if(n == 0) {
      this->n = field.getNbComponent() - stride;
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  
  virtual void registerToDumper(const std::string & id, iohelper::Dumper & dumper) {
    dumper.addNodeDataField(id, *this);
  }

  inline iterator begin() {
    return iterator(field.storage(), field.getNbComponent(), 
		    n, stride, filter->storage());
  }

  inline iterator end() {
    return iterator(field.storage(), field.getNbComponent(), 
		    n, stride, filter->storage()+filter->getSize());
  }

  bool isHomogeneous() {
    return true;
  }
  void checkHomogeneity() { this->homogeneous = true; }

  virtual UInt getDim() {
    if(this->padding) return this->padding;
    else              return n;
  }

  void setPadding(UInt padding){this->padding = padding;}

  UInt size() {
    return filter->getSize();
  }

  iohelper::DataType getDataType() {
    return iohelper::getDataType<T>();
  }


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  
private:
  const Container & field;
  UInt n, stride;
  const Filter * filter;

  UInt padding;
};

 
__END_AKANTU_DUMPER__
__END_AKANTU__
/* -------------------------------------------------------------------------- */
#endif /* __AKANTU_DUMPER_NODAL_FIELD_HH__ */
