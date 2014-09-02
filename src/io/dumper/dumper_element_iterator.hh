#ifndef __AKANTU_DUMPER_ELEMENT_ITERATOR_HH__
#define __AKANTU_DUMPER_ELEMENT_ITERATOR_HH__
/* -------------------------------------------------------------------------- */
#include "element.hh"
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */

template<class types, template <class> class final_iterator>
class element_iterator : 
  public iohelper::iterator< 
                            typename types::data_type, 
                            final_iterator<types>,
                            typename types::return_type>

{

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
  
public:

  typedef typename types::it_type      it_type;
  typedef typename types::field_type field_type; 
  typedef typename types::array_type array_type; 
  typedef typename types::array_iterator array_iterator; 
  typedef final_iterator<types> iterator;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  element_iterator(const field_type & field,
		   const typename field_type::type_iterator & t_it,
                   const typename field_type::type_iterator & t_it_end,
                   const array_iterator & array_it,
                   const array_iterator & array_it_end,
                   const GhostType ghost_type = _not_ghost) 
    : field(field),
      tit(t_it),
      tit_end(t_it_end),
      array_it(array_it),
      array_it_end(array_it_end),
      ghost_type(ghost_type) {
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  
public:

  bool operator!=(const iterator & it) const {
    return (ghost_type != it.ghost_type) 
      || (tit != it.tit || (array_it != it.array_it));
  }

  iterator & operator++() {
    ++array_it;
    while(array_it == array_it_end && tit != tit_end) {
      ++tit;
      if(tit != tit_end) {

	const array_type & vect = field(*tit, ghost_type);
	UInt _nb_data_per_elem = getNbDataPerElem(*tit);
	UInt nb_component = vect.getNbComponent();
	UInt size = (vect.getSize() * nb_component) / _nb_data_per_elem;

	array_it       = vect.begin_reinterpret(_nb_data_per_elem,size);
        array_it_end   = vect.end_reinterpret  (_nb_data_per_elem,size);
      }
    }
    return *(static_cast<iterator *>(this));
  };

  ElementType getType() { return *tit; }

  Element getCurrentElement(){
    return Element(*tit,array_it.getCurrentIndex());
  }

  UInt getNbDataPerElem(const ElementType & type) const { 
    return nb_data_per_elem(type,ghost_type);
  }
  
  void setNbDataPerElem(const ElementTypeMap<UInt> & nb_data){
    this->nb_data_per_elem = nb_data;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  
protected:

  /// the field to iterate on
  const field_type & field;
  /// field iterator
  typename field_type::type_iterator tit;
  /// field iterator end
  typename field_type::type_iterator tit_end;
  /// array iterator
  array_iterator array_it;
  /// internal iterator end
  array_iterator array_it_end;
  /// ghost type identification
  const GhostType ghost_type;
  /// number of data per element
  ElementTypeMap<UInt> nb_data_per_elem;

};

/* -------------------------------------------------------------------------- */
template<typename types>
class elemental_field_iterator 
  : public element_iterator<types, elemental_field_iterator> {

public:

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  

  typedef element_iterator<types,elemental_field_iterator> parent;
  typedef typename types::it_type     it_type;
  typedef typename types::return_type return_type;
  typedef typename types::field_type  field_type;
  typedef typename types::array_iterator array_iterator;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  elemental_field_iterator(const field_type & field,
			   const typename field_type::type_iterator & t_it,
                           const typename field_type::type_iterator & t_it_end,
                           const array_iterator & array_it,
			   const array_iterator & array_it_end,
			   const GhostType ghost_type = _not_ghost) :
    parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) { }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  
  return_type operator*(){
    return *this->array_it;
  }

private:

};

/* -------------------------------------------------------------------------- */
__END_AKANTU_DUMPER__
__END_AKANTU__
/* -------------------------------------------------------------------------- */



#endif /* __AKANTU_DUMPER_ELEMENT_ITERATOR_HH__ */
