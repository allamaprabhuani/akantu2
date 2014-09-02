#include "dumper_generic_elemental_field.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */

template <class types>
class filtered_connectivity_field_iterator 
  : public element_iterator<types,filtered_connectivity_field_iterator> {

public:

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
  

  typedef element_iterator<types,filtered_connectivity_field_iterator> parent;
  typedef typename types::return_type return_type;
  typedef typename types::field_type field_type;
  typedef typename types::array_iterator array_iterator;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  filtered_connectivity_field_iterator(const field_type & field,
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
    const Vector<UInt> & old_connect = *this->array_it;
    Vector<UInt> new_connect(old_connect.size());
    Array<UInt>::const_iterator<UInt> nodes_begin = nodal_filter->begin();
    Array<UInt>::const_iterator<UInt> nodes_end = nodal_filter->end();
    for(UInt i(0); i<old_connect.size(); ++i) {
      Array<UInt>::const_iterator<UInt> new_id =
        std::find(nodes_begin, nodes_end, old_connect(i));
      if(new_id == nodes_end) AKANTU_EXCEPTION("Node not found in the filter!");
      new_connect(i) = new_id - nodes_begin;
    }
    return new_connect;
  }

  void setNodalFilter(const Array<UInt> & new_nodal_filter) {
    nodal_filter = &new_nodal_filter;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  
private:
  const Array<UInt> * nodal_filter;
};

/* -------------------------------------------------------------------------- */

class FilteredConnectivityField : 
  public GenericElementalField<SingleType<UInt,Vector,true>,
 			       filtered_connectivity_field_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  

public:

  typedef SingleType<UInt,Vector,true> types;
  typedef filtered_connectivity_field_iterator<types> iterator;
  typedef types::field_type field_type;
  typedef GenericElementalField<types,filtered_connectivity_field_iterator> parent;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  FilteredConnectivityField(const field_type & field,
 			    const Array<UInt> & nodal_filter,
			    UInt spatial_dimension = _all_dimensions,
			    GhostType ghost_type = _not_ghost,
			    ElementKind element_kind = _ek_not_defined) :
    parent(field, spatial_dimension, ghost_type, element_kind),
    nodal_filter(nodal_filter) { }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  
  iterator begin() {
    iterator it = parent::begin();
    it.setNodalFilter(nodal_filter);
    return it;
  }

  iterator end() {
    iterator it = parent::end();
    it.setNodalFilter(nodal_filter);
    return it;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */


private:
  const Array<UInt> & nodal_filter;
};


/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
__END_AKANTU__

/* -------------------------------------------------------------------------- */
