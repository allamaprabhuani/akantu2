
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */


template<class types>
class element_partition_field_iterator 
  : public element_iterator<types,element_partition_field_iterator> {

public:
  
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
  
  typedef element_iterator<types, element_partition_field_iterator> parent;
  typedef typename types::return_type return_type;
  typedef typename types::array_iterator array_iterator;
  typedef typename types::field_type field_type;


public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  element_partition_field_iterator(const field_type & field,
				   const typename field_type::type_iterator & t_it,
                                   const typename field_type::type_iterator & t_it_end,
                                   const array_iterator & array_it,
				   const array_iterator & array_it_end,
                                   const GhostType ghost_type = _not_ghost) :
    parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) {
    prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  
  return_type operator*() {
    return return_type(1, prank);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  
protected:
  UInt prank;
};

/* -------------------------------------------------------------------------- */

template<bool filtered = false>
class ElementPartitionField : 
  public GenericElementalField<SingleType<UInt,Vector,filtered>,
 			       element_partition_field_iterator> {
public:

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  

  typedef SingleType<UInt,Vector,filtered> types;
  typedef element_partition_field_iterator<types> iterator;
  typedef GenericElementalField<types,element_partition_field_iterator> parent;
  typedef typename types::field_type field_type;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  ElementPartitionField(const field_type & field,
			UInt spatial_dimension = _all_dimensions,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    parent(field, spatial_dimension, ghost_type, element_kind) {
    this->homogeneous = true;
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  UInt getDim() { return 1; }
};

/* -------------------------------------------------------------------------- */
__END_AKANTU_DUMPER__
__END_AKANTU__
