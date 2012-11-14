/**
 * @file   dumper_iohelper_tmpl_elemental_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 26 21:52:40 2012
 *
 * @brief  description of elemental fields
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

__END_AKANTU__
#include "static_communicator.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename i_type, typename d_type,
	 template<typename> class ret_type, class daughter>
class DumperIOHelper::element_iterator : public iohelper::iterator< d_type, daughter,
								    ret_type<d_type> > {
public:
  typedef i_type              it_type;
  typedef d_type              data_type;
  typedef ret_type<data_type> return_type;
  typedef ByElementTypeVector<it_type> field_type;
  typedef typename Vector<it_type>::template const_iterator< ret_type<it_type> > internal_iterator;
public:
  element_iterator(const field_type & field,
		   UInt n,
		   const typename field_type::type_iterator & t_it,
		   const typename field_type::type_iterator & t_it_end,
		   const internal_iterator & it,
		   ElementType element_type,
		   const GhostType ghost_type = _not_ghost) : field(field),
							      n(n),
							      padding_n(0),
							      padding_m(0),
							      itn(0), itm(0),
							      tit(t_it),
							      tit_end(t_it_end),
							      vit(it),
							      ghost_type(ghost_type) {
    UInt nb_data = getNbDataPerElem(element_type);
    const Vector<it_type> & vect = field(element_type, ghost_type);
    UInt size = vect.getSize() / nb_data;
    UInt ln = n;
    if(n == 0) ln = vect.getNbComponent();
    vit_end   = iterator_helper<it_type, ret_type>::end(vect, ln, vect.getNbComponent() / ln * nb_data, size);
  }

public:
  bool operator!=(const daughter & it) const {
    return (ghost_type != ghost_type) || (tit != it.tit || (vit != it.vit));
  }

  daughter & operator++() {
    ++vit;
    if(vit == vit_end) {
      ++tit;
      if(tit != tit_end) {
	UInt nb_data = getNbDataPerElem(*tit);
	const Vector<it_type> & vect = field(*tit, ghost_type);
	UInt size = vect.getSize() / nb_data;
	UInt ln = n, lm;
	if(n == 0) ln = vect.getNbComponent();
	if(itn != 0) ln = itn;
	if(itm != 0) lm = itm;
	else lm = vect.getNbComponent() / ln * nb_data;
	vit     = iterator_helper<it_type, ret_type>::begin(vect, ln, lm, size);
	vit_end = iterator_helper<it_type, ret_type>::end  (vect, ln, lm, size);
      }
    }
    return *(static_cast<daughter *>(this));
  };

  ElementType getType() { return *tit; }

  void setPadding(UInt n, UInt m) {
    padding_n = n;
    padding_m = m;
  }

  // small trick for the material iterator
  void setItSize(UInt n, UInt m) { itn = n; itm = m; }

protected:
  virtual UInt getNbDataPerElem(__attribute__((unused)) const ElementType & type) { return 1; }

protected:
  const field_type & field;
  UInt n, padding_n, padding_m;
  UInt itn, itm;
  typename field_type::type_iterator tit;
  typename field_type::type_iterator tit_end;
  internal_iterator vit;
  internal_iterator vit_end;
  const GhostType ghost_type;
};


/* -------------------------------------------------------------------------- */
template<typename T, template<class> class ret_type>
class DumperIOHelper::elemental_field_iterator : public element_iterator< T, T, ret_type,
									  elemental_field_iterator<T, ret_type> > {
public:
  typedef element_iterator< T, T, ret_type,
			    elemental_field_iterator<T, ret_type> > parent;
  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
public:
  elemental_field_iterator(const field_type & field,
			   UInt n,
			   const typename field_type::type_iterator & t_it,
			   const typename field_type::type_iterator & t_it_end,
			   const internal_iterator & it,
			   ElementType element_type,
			   const GhostType ghost_type = _not_ghost) :
    parent(field, n, t_it, t_it_end, it, element_type, ghost_type) { }

  return_type operator*(){
    return padding_helper.pad(*this->vit, this->padding_m, this->padding_n, this->getNbDataPerElem(*this->tit));
  }
private:
  PaddingHelper<T, ret_type> padding_helper;
};

/* -------------------------------------------------------------------------- */
class DumperIOHelper::element_type_field_iterator : public element_iterator<UInt, iohelper::ElemType,
									    types::Vector,
									    element_type_field_iterator> {
public:
  typedef element_iterator<UInt, iohelper::ElemType,
			   types::Vector, element_type_field_iterator> parent;

  typedef parent::it_type     it_type;
  typedef parent::data_type   data_type;
  typedef parent::return_type return_type;
  typedef parent::field_type  field_type;
  typedef parent::internal_iterator internal_iterator;
public:
  element_type_field_iterator(const field_type & field,
			      __attribute__((unused)) UInt n,
			      const field_type::type_iterator & t_it,
			      const field_type::type_iterator & t_it_end,
			      const internal_iterator & it,
			      ElementType element_type,
			      const GhostType ghost_type = _not_ghost) :
    parent(field, 0, t_it, t_it_end, it, element_type, ghost_type) { }

  return_type operator*() {
    data_type type = DumperIOHelper::getIOHelperType(*tit);
    return return_type(1, type);
  }
};


/* -------------------------------------------------------------------------- */
class DumperIOHelper::element_partition_field_iterator : public element_iterator<UInt, UInt,
										 types::Vector,
										 element_partition_field_iterator> {
public:
  typedef element_iterator<UInt, UInt,
			   types::Vector, element_partition_field_iterator> parent;

  typedef parent::it_type     it_type;
  typedef parent::data_type   data_type;
  typedef parent::return_type return_type;
  typedef parent::field_type  field_type;
  typedef parent::internal_iterator internal_iterator;
public:
  element_partition_field_iterator(const field_type & field,
				   UInt n,
				   const field_type::type_iterator & t_it,
				   const field_type::type_iterator & t_it_end,
				   const internal_iterator & it,
				   ElementType element_type,
				   const GhostType ghost_type = _not_ghost) :
    parent(field, n, t_it, t_it_end, it, element_type, ghost_type) {
    prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  }

  return_type operator*() {
    return return_type(1, prank);
  }
protected:
  UInt prank;
};

/* -------------------------------------------------------------------------- */
/* Fields type description                                                    */
/* -------------------------------------------------------------------------- */
template<typename T, class iterator_type, template<class> class ret_type>
class DumperIOHelper::GenericElementalField : public Field {
public:
  GenericElementalField(const ByElementTypeVector<T> & field,
			UInt spatial_dimension = 0,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    field(field), spatial_dimension(spatial_dimension),
    ghost_type(ghost_type), element_kind(element_kind), n(0), itn(0) {
    UInt nb_component;
    homogeneous = checkHomogeneity(field, nb_component, nb_total_element);
    if(homogeneous && n == 0) n = nb_component;
  }

  GenericElementalField(const ByElementTypeVector<T> & field,
			UInt n,
			UInt spatial_dimension = 0,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    field(field), spatial_dimension(spatial_dimension),
    ghost_type(ghost_type), element_kind(element_kind), n(n), itn(0) {
    UInt nb_component;
    homogeneous = checkHomogeneity(field, nb_component, nb_total_element);
  }

  typedef iterator_type iterator;
  typedef typename iterator::internal_iterator internal_iterator;
  typedef typename ByElementTypeVector<T>::type_iterator type_iterator;
  /* ------------------------------------------------------------------------ */
  virtual iterator begin() {
    type_iterator tit = field.firstType(spatial_dimension, ghost_type, element_kind);
    type_iterator end = field.lastType(spatial_dimension, ghost_type, element_kind);

    const Vector<T> & vect = field(*tit, ghost_type);
    UInt nb_data = getNbDataPerElem(*tit);
    UInt nb_component = vect.getNbComponent();
    UInt size = vect.getSize() / nb_data;
    UInt ln = nb_component;
    if(itn != 0) ln = itn;

    internal_iterator it = iterator_helper<T, ret_type>::begin(vect, ln, nb_component / ln * nb_data, size);
    iterator rit = iterator(field, n, tit, end, it, *tit, ghost_type);
    rit.setPadding(padding_n, padding_m);
    return rit;
  }

  virtual iterator end  () {
    type_iterator tit = field.firstType(spatial_dimension, ghost_type, element_kind);
    type_iterator end = field.lastType(spatial_dimension, ghost_type, element_kind);

    ElementType type = *tit;
    for (; tit != end; ++tit) type = *tit;

    const Vector<T> & vect = field(type, ghost_type);
    UInt nb_data = getNbDataPerElem(type);
    UInt nb_component = vect.getNbComponent();
    UInt size = vect.getSize() / nb_data;
    UInt ln = nb_component;
    if(itn != 0 ) ln = itn;
    internal_iterator it = iterator_helper<T, ret_type>::end(vect, ln, nb_component / ln * nb_data, size);
    iterator rit =  iterator(field, n, end, end, it, type, ghost_type);
    rit.setPadding(padding_n, padding_m);
    return rit;

  }

  virtual void registerToDumper(const std::string & id, iohelper::Dumper & dupmer) {
    dupmer.addElemDataField(id, *this);
  };

  bool isHomogeneous() { return homogeneous; }
  virtual UInt getDim() {
    if(padding_n && padding_m)
      return padding_m*padding_n;
    else return n;
  }
  UInt size() { return nb_total_element; }

protected:
  virtual UInt getNbDataPerElem(__attribute__((unused)) const ElementType & type) { return 1; }

  virtual bool checkHomogeneity(const ByElementTypeVector<T> & field_to_check,
				UInt & nb_comp, UInt & nb_elem) {
    type_iterator tit = field_to_check.firstType(spatial_dimension, ghost_type, element_kind);
    type_iterator end = field_to_check.lastType(spatial_dimension, ghost_type, element_kind);
    nb_elem = 0;
    nb_comp = 0;
    UInt nb_data = getNbDataPerElem(*tit);

    bool homogen = true;
    if(tit != end) {
      nb_comp = field_to_check(*tit, ghost_type).getNbComponent() * nb_data;
      for(;tit != end; ++tit) {
	const Vector<T> & vect = field_to_check(*tit, ghost_type);
	UInt nb_data_cur = getNbDataPerElem(*tit);
	UInt nb_element = vect.getSize() / nb_data_cur;
	UInt nb_comp_cur = vect.getNbComponent() * nb_data_cur;
	if(homogen && nb_comp != nb_comp_cur) homogen = false;
	nb_elem += nb_element;
      }

      if(!homogen) nb_comp = 0;
    }

    return homogen;
  }

protected:
  const ByElementTypeVector<T> & field;
  UInt nb_total_element;
  UInt spatial_dimension;
  GhostType ghost_type;
  ElementKind element_kind;
  bool homogeneous;
  UInt n, itn;
};

/* -------------------------------------------------------------------------- */
template<typename T, template<typename> class ret_type>
class DumperIOHelper::ElementalField : public GenericElementalField<T, elemental_field_iterator<T, ret_type>, ret_type> {
public:
  typedef elemental_field_iterator<T, ret_type> iterator;

  ElementalField(const ByElementTypeVector<T> & field,
		 UInt spatial_dimension = 0,
		 GhostType ghost_type = _not_ghost,
		 ElementKind element_kind = _ek_not_defined) :
    GenericElementalField<T, iterator, ret_type>(field, 0, spatial_dimension,
						 ghost_type, element_kind) { }
};

template<typename T>
class DumperIOHelper::ElementalField<T, types::Matrix> : public GenericElementalField<T, elemental_field_iterator<T, types::Matrix>, types::Matrix> {
public:
  typedef elemental_field_iterator<T, types::Matrix> iterator;

  ElementalField(const ByElementTypeVector<T> & field,
		 UInt n,
		 UInt spatial_dimension = 0,
		 GhostType ghost_type = _not_ghost,
		 ElementKind element_kind = _ek_not_defined) :
    GenericElementalField<T, iterator, types::Matrix>(field, n, spatial_dimension,
						      ghost_type, element_kind) { }
};

/* -------------------------------------------------------------------------- */
class DumperIOHelper::ElementTypeField : public GenericElementalField<UInt,
								      element_type_field_iterator,
								      types::Vector> {
public:
  typedef element_type_field_iterator iterator;
private:
  typedef GenericElementalField<UInt, iterator, types::Vector> parent;
public:
  /* ------------------------------------------------------------------------ */
  ElementTypeField(const Mesh & mesh,
		   UInt spatial_dimension = 0,
		   GhostType ghost_type = _not_ghost,
		   ElementKind element_kind = _ek_not_defined) :
    parent(mesh.getConnectivities(), 0, spatial_dimension, ghost_type, element_kind) {
    homogeneous = true;
  }

  virtual void registerToDumper(__attribute__((unused)) const std::string & id,
			   iohelper::Dumper & dupmer) {
    dupmer.addElemDataField("element_type", *this);
  }

  UInt getDim() { return 1; }
};


/* -------------------------------------------------------------------------- */
class DumperIOHelper::ElementPartitionField : public GenericElementalField<UInt,
									   element_partition_field_iterator,
									   types::Vector> {
public:
  typedef element_partition_field_iterator iterator;
private:
  typedef GenericElementalField<UInt, iterator, types::Vector> parent;
public:
  /* ------------------------------------------------------------------------ */
  ElementPartitionField(const Mesh & mesh,
			UInt spatial_dimension = 0,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    parent(mesh.getConnectivities(), spatial_dimension, ghost_type, element_kind) {
    homogeneous = true;
  }

  UInt getDim() { return 1; }
};

