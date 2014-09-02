#ifndef __AKANTU_DUMPER_TYPE_TRAITS_HH__
#define __AKANTU_DUMPER_TYPE_TRAITS_HH__
/* -------------------------------------------------------------------------- */
#include "element_type_map.hh"
#include "element_type_map_filter.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */

template <class data, class ret, class field>
struct TypeTraits {

  //! the stored data (real, int, uint, ...)
  typedef data  data_type;
  //! the type returned by the operator *
  typedef ret   return_type;
  //! the field type (ElementTypeMap or ElementTypeMapFilter)
  typedef field field_type;
  //! the type over which we iterate
  typedef typename field_type::type it_type;
  //! the type of array (Array<T> or ArrayFilter<T>)
  typedef typename field_type::array_type array_type;
  //! the iterator over the array
  typedef typename array_type::const_vector_iterator array_iterator;
};

/* -------------------------------------------------------------------------- */

template <class T, template <class> class ret, bool filtered>
struct SingleType
  : public TypeTraits<T,ret<T>,ElementTypeMapArray<T> >{
  
};

/* -------------------------------------------------------------------------- */
template <class T, template <class> class ret>
struct SingleType<T,ret,true> : 
  public TypeTraits<T,ret<T>, ElementTypeMapArrayFilter<T> >{
  
};
/* -------------------------------------------------------------------------- */
template <class it_type, class data_type, template <class> class ret, 
	  bool filtered>
struct DualType
  : public TypeTraits<data_type,ret<data_type>, ElementTypeMapArray<it_type> >{
};

/* -------------------------------------------------------------------------- */
template <class it_type, class data_type,template <class> class ret>
struct DualType<it_type,data_type,ret,true> : 
  public TypeTraits<data_type,ret<data_type>, ElementTypeMapArrayFilter<it_type> >{
  
};
/* -------------------------------------------------------------------------- */
__END_AKANTU_DUMPER__
__END_AKANTU__

/* -------------------------------------------------------------------------- */


#endif /* __AKANTU_DUMPER_TYPE_TRAITS_HH__ */
