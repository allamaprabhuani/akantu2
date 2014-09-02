/* -------------------------------------------------------------------------- */

template<class types, template <class> class iterator>
void GenericElementalField<types,iterator>
::checkHomogeneity() {
  

  typedef typename field_type::type_iterator field_type_iterator;
  typedef typename field_type::array_type array_type;

  field_type_iterator tit;
  field_type_iterator end;
  
  tit = field.firstType(spatial_dimension, ghost_type, element_kind);
  end = field.lastType(spatial_dimension, ghost_type, element_kind);
  
  this->nb_total_element = 0;
  UInt nb_comp = 0;
  
  bool homogen = true;
  if(tit != end) {
    nb_comp = this->field(*tit, ghost_type).getNbComponent();
    for(;tit != end; ++tit) {
      const array_type & vect = this->field(*tit, ghost_type);
      UInt nb_element = vect.getSize();
      UInt nb_comp_cur = vect.getNbComponent();
      if(homogen && nb_comp != nb_comp_cur) homogen = false;
      this->nb_total_element += nb_element;

      //      this->nb_data_per_elem(*tit,this->ghost_type) = nb_comp_cur;
    }
    
    if(!homogen) nb_comp = 0;
  }
  
  this->homogeneous = homogen;
}

/* -------------------------------------------------------------------------- */


