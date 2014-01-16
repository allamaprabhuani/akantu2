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
#include "mesh.hh"
#include "static_communicator.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename i_type, typename d_type,
         template<typename> class ret_type, class daughter>
class DumperIOHelper::element_iterator<i_type, d_type,
                                       ret_type, daughter,
                                       false> : public iohelper::iterator< d_type,
                                                                           daughter,
                                                                           ret_type<d_type> > {
public:
  typedef i_type              it_type;
  typedef d_type              data_type;
  typedef ret_type<data_type> return_type;
  typedef ByElementTypeArray<it_type> field_type;
  typedef typename Array<it_type>::template const_iterator< ret_type<it_type> > internal_iterator;
public:
  element_iterator(const field_type & field,
                   UInt n,
                   const typename field_type::type_iterator & t_it,
                   const typename field_type::type_iterator & t_it_end,
                   const internal_iterator & it,
                   ElementType element_type,
                   const GhostType ghost_type = _not_ghost,
                   const ByElementTypeArray<UInt> * filter = NULL,
                   UInt * fit = NULL) : field(field),
                                        n(n),
                                        padding_n(0),
                                        padding_m(0),
                                        itn(0), itm(0),
                                        tit(t_it),
                                        tit_end(t_it_end),
                                        vit(it),
                                        ghost_type(ghost_type),
    filter(filter),
    fit(fit) {
    UInt nb_data = getNbDataPerElem(element_type);
    const Array<it_type> & vect = field(element_type, ghost_type);
    UInt size = vect.getSize() / nb_data;
    UInt ln = n;
    if(n == 0) ln = vect.getNbComponent();
    vit_end   = iterator_helper<it_type, ret_type>::end(vect,
                                                        ln,
                                                        vect.getNbComponent() / ln * nb_data, size);
  }

public:
  bool operator!=(const daughter & it) const {
    return (ghost_type != it.ghost_type) || (tit != it.tit || (vit != it.vit));
  }

  daughter & operator++() {
    ++vit;
    while(vit == vit_end && tit != tit_end) {
      ++tit;
      if(tit != tit_end) {
        UInt nb_data = getNbDataPerElem(*tit);
        const Array<it_type> & vect = field(*tit, ghost_type);
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
  const ByElementTypeArray<UInt> * filter;
  UInt * fit;
};


/* -------------------------------------------------------------------------- */
template<typename i_type, typename d_type,
         template<typename> class ret_type,
         class daughter>
class DumperIOHelper::element_iterator<i_type, d_type, ret_type, daughter, true> : public iohelper::iterator< d_type, daughter, ret_type<d_type> > {
public:
  typedef i_type              it_type;
  typedef d_type              data_type;
  typedef ret_type<data_type> return_type;
  typedef ByElementTypeArray<it_type> field_type;
  typedef typename Array<it_type>::template const_iterator< ret_type<it_type> > internal_iterator;
public:
  element_iterator(const field_type & field,
                   UInt n,
                   const typename field_type::type_iterator & t_it,
                   const typename field_type::type_iterator & t_it_end,
                   const internal_iterator & it,
                   ElementType element_type,
                   const GhostType ghost_type = _not_ghost,
                   const ByElementTypeArray<UInt> * filter = NULL,
                   UInt * fit = NULL) : field(field),
                                        n(n),
                                        padding_n(0),
                                        padding_m(0),
                                        itn(0), itm(0),
                                        tit(t_it),
                                        tit_end(t_it_end),
                                        vit(it),
                                        ghost_type(ghost_type),
    filter(filter),
    fit(fit) {
    UInt nb_data = getNbDataPerElem(element_type);
    const Array<it_type> & vect = field(element_type, ghost_type);
    UInt size = vect.getSize() / nb_data;
    UInt ln = n;
    if(n == 0) ln = vect.getNbComponent();
    UInt lm = vect.getNbComponent() / ln * nb_data;
    const Array<UInt> & filter_array = (*filter)(element_type, ghost_type);

    this->fit_end = filter_array.storage() + filter_array.getSize();
    this->vit_begin = iterator_helper<it_type, ret_type>::begin(vect, ln, lm, size);

    searchNonEmptyArray();
    if(this->fit != this->fit_end) this->vit = this->vit_begin + *this->fit;
  }

public:
  bool operator!=(const daughter & it) const {
    return (ghost_type != it.ghost_type) || (tit != it.tit || (fit != it.fit));
  }

  daughter & operator++() {
    ++fit;
    searchNonEmptyArray();

    if(fit != fit_end) vit = vit_begin + *fit;
    return *(static_cast<daughter *>(this));
  };

  ElementType getType() { return *tit; }

  void setPadding(UInt n, UInt m) {
    padding_n = n;
    padding_m = m;
  }

  // small trick for the material iterator
  void setItSize(UInt n, UInt m) {
    itn = n;
    itm = m;
    if(tit != tit_end) {
      const Array<it_type> & vect = field(*tit, ghost_type);
      UInt nb_data = getNbDataPerElem(*tit);
      UInt size = vect.getSize() / nb_data;
      vit_begin = iterator_helper<it_type, ret_type>::begin(vect, itn, itm, size);
      if(fit != fit_end) vit = vit_begin + *fit;
    }
  }

protected:
  void searchNonEmptyArray() {
    while(fit == fit_end && tit != tit_end) {
      ++tit;
      if(tit != tit_end) {
        UInt nb_data = getNbDataPerElem(*tit);
        const Array<it_type> & vect = field(*tit, ghost_type);
        UInt size = vect.getSize() / nb_data;
        UInt ln = n, lm;
        if(n == 0) ln = vect.getNbComponent();
        if(itn != 0) ln = itn;
        if(itm != 0) lm = itm;
        else lm = vect.getNbComponent() / ln * nb_data;
        const Array<UInt> & f = (*filter)(*tit, ghost_type);
        fit     = f.storage();
        fit_end = fit + f.getSize();
        vit_begin = iterator_helper<it_type, ret_type>::begin(vect, ln, lm, size);
      }
    }
  }

protected:
  virtual UInt getNbDataPerElem(__attribute__((unused)) const ElementType & type) { return 1; }

protected:
  const field_type & field;
  UInt n, padding_n, padding_m;
  UInt itn, itm;
  typename field_type::type_iterator tit;
  typename field_type::type_iterator tit_end;
  internal_iterator vit_begin;
  internal_iterator vit;
  const GhostType ghost_type;
  const ByElementTypeArray<UInt> * filter;
  UInt * fit, * fit_end;
};


/* -------------------------------------------------------------------------- */
template<typename T, template<class> class ret_type, bool filtered>
class DumperIOHelper::elemental_field_iterator : public element_iterator< T, T, ret_type,
                                                                          elemental_field_iterator<T, ret_type, filtered>, filtered > {
public:
  typedef element_iterator< T, T, ret_type,
                            elemental_field_iterator<T, ret_type, filtered>, filtered > parent;
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
                           const GhostType ghost_type = _not_ghost,
                           const ByElementTypeArray<UInt> * filter = NULL,
                           UInt * fit = NULL) :
    parent(field, n, t_it, t_it_end, it, element_type, ghost_type, filter, fit) { }

  return_type operator*(){
    return padding_helper.pad(*this->vit, this->padding_m, this->padding_n, this->getNbDataPerElem(*this->tit));
  }
private:
  PaddingHelper<T, ret_type> padding_helper;
};

/* -------------------------------------------------------------------------- */
class DumperIOHelper::filtered_connectivity_field_iterator : public element_iterator<UInt, UInt, Vector,
                                                                                     filtered_connectivity_field_iterator, true> {
public:
  typedef element_iterator<UInt, UInt, Vector, filtered_connectivity_field_iterator, true> parent;
  typedef parent::it_type     it_type;
  typedef parent::data_type   data_type;
  typedef parent::return_type return_type;
  typedef parent::field_type  field_type;
  typedef parent::internal_iterator internal_iterator;
public:
  filtered_connectivity_field_iterator(const field_type & field,
                                       UInt n,
                                       const field_type::type_iterator & t_it,
                                       const field_type::type_iterator & t_it_end,
                                       const internal_iterator & it,
                                       ElementType element_type,
                                       const GhostType ghost_type = _not_ghost,
                                       const ByElementTypeArray<UInt> * filter = NULL,
                                       UInt * fit = NULL) :
  parent(field, n, t_it, t_it_end, it, element_type, ghost_type, filter, fit) { }

  return_type operator*(){
    const Vector<UInt> & old_connect = *this->vit;
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

private:
  const Array<UInt> * nodal_filter;
};

/* -------------------------------------------------------------------------- */
template<bool filtered>
class DumperIOHelper::element_type_field_iterator : public element_iterator<UInt, iohelper::ElemType,
                                                                            Vector,
                                                                            element_type_field_iterator<filtered>,
                                                                            filtered> {
public:
  typedef element_iterator<UInt, iohelper::ElemType,
                           Vector, element_type_field_iterator<filtered>, filtered> parent;

  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
public:
  element_type_field_iterator(const field_type & field,
                              __attribute__((unused)) UInt n,
                              const typename field_type::type_iterator & t_it,
                              const typename field_type::type_iterator & t_it_end,
                              const internal_iterator & it,
                              ElementType element_type,
                              const GhostType ghost_type = _not_ghost,
                              const ByElementTypeArray<UInt> * filter = NULL,
                              UInt * fit = NULL) :
    parent(field, 0, t_it, t_it_end, it, element_type, ghost_type, filter, fit) { }

  return_type operator*() {
    data_type type = getIOHelperType(*this->tit);
    return return_type(1, type);
  }
};


/* -------------------------------------------------------------------------- */
template<bool filtered>
class DumperIOHelper::element_partition_field_iterator : public element_iterator<UInt, UInt,
                                                                                 Vector,
                                                                                 element_partition_field_iterator<filtered>, filtered> {
public:
  typedef element_iterator<UInt, UInt,
                           Vector, element_partition_field_iterator<filtered>, filtered > parent;

  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
public:
  element_partition_field_iterator(const field_type & field,
                                   __attribute__((unused)) UInt n,
                                   const typename field_type::type_iterator & t_it,
                                   const typename field_type::type_iterator & t_it_end,
                                   const internal_iterator & it,
                                   ElementType element_type,
                                   const GhostType ghost_type = _not_ghost,
                                   const ByElementTypeArray<UInt> * filter = NULL,
                                   UInt * fit = NULL) :
    parent(field, 0, t_it, t_it_end, it, element_type, ghost_type, filter, fit) {
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
template<typename T, class iterator_type, template<class> class ret_type, bool filtered>
class DumperIOHelper::GenericElementalField : public Field {
public:
  GenericElementalField(const ByElementTypeArray<T> & field,
                        UInt spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_not_defined,
                        const ByElementTypeArray<UInt> * filter = NULL) :
    field(field), spatial_dimension(spatial_dimension),
    ghost_type(ghost_type), element_kind(element_kind), n(0), itn(0), filter(filter) {
    AKANTU_DEBUG_ASSERT(filtered == (filter != NULL) , "Filter problem!");
    UInt nb_component;
    homogeneous = checkHomogeneity(field, nb_component, nb_total_element);
    out_n = nb_component;
  }

  GenericElementalField(const ByElementTypeArray<T> & field,
                        UInt n,
                        UInt spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_not_defined,
                        const ByElementTypeArray<UInt> * filter = NULL) :
    field(field), spatial_dimension(spatial_dimension),
    ghost_type(ghost_type), element_kind(element_kind), n(n), itn(0), filter(filter) {
    AKANTU_DEBUG_ASSERT(filtered == (filter != NULL) , "Filter problem!");
    UInt nb_component;
    homogeneous = checkHomogeneity(field, nb_component, nb_total_element);
    out_n = nb_component;
  }

  typedef iterator_type iterator;
  typedef typename iterator::internal_iterator internal_iterator;
  typedef typename ByElementTypeArray<T>::type_iterator type_iterator;
  /* ------------------------------------------------------------------------ */
  virtual iterator begin() {
    type_iterator tit;
    type_iterator end;

    /// type iterators on the elemental field
    tit = field.firstType(spatial_dimension, ghost_type, element_kind);
    end = field.lastType(spatial_dimension, ghost_type, element_kind);

    /// skip all types without data
    ElementType type = *tit;
    for (;tit != end && field(*tit, ghost_type).getSize() == 0; ++tit)
      type = *tit;

    /// getting information for the field of the given type
    const Array<T> & vect = field(type, ghost_type);
    UInt * fit = NULL;
    if(filtered) {
      const Array<UInt> & typed_filter = (*filter)(*tit, ghost_type);
      fit = typed_filter.storage();
    }
    UInt nb_data = getNbDataPerElem(*tit);
    UInt nb_component = vect.getNbComponent();
    UInt size = vect.getSize() / nb_data;
    UInt ln = nb_component;

    /// if nb_components to be dumped is explicitly defined, use it
    if(itn != 0) ln = itn;

    /// define element-wise iterator
    internal_iterator it = iterator_helper<T, ret_type>::begin(vect, ln, nb_component / ln * nb_data, size);

    /// define data iterator
    iterator rit = iterator(field, n, tit, end, it, type, ghost_type, filter, fit);

    /// set user defined padding
    rit.setPadding(padding_n, padding_m);

    return rit;
  }

  virtual iterator end  () {
    type_iterator tit;
    type_iterator end;

    tit = field.firstType(spatial_dimension, ghost_type, element_kind);
    end = field.lastType(spatial_dimension, ghost_type, element_kind);


    ElementType type = *tit;
    for (; tit != end; ++tit) type = *tit;

    const Array<T> & vect = field(type, ghost_type);
    UInt * fit = NULL;
    if(filtered) {
      const Array<UInt> & typed_filter = (*filter)(type, ghost_type);
      fit = typed_filter.storage()+typed_filter.getSize();
    }
    UInt nb_data = getNbDataPerElem(type);
    UInt nb_component = vect.getNbComponent();
    UInt size = vect.getSize() / nb_data;
    UInt ln = nb_component;
    if(itn != 0 ) ln = itn;
    internal_iterator it = iterator_helper<T, ret_type>::end(vect, ln, nb_component / ln * nb_data, size);
    iterator rit =  iterator(field, n, end, end, it, type, ghost_type, filter, fit);
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
    else if(n == 0) return out_n;
    else return n;
  }
  UInt size() {
    UInt nb_component;
    checkHomogeneity(field, nb_component, nb_total_element);
    return nb_total_element;
  }

  iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }

protected:
  virtual UInt getNbDataPerElem(__attribute__((unused)) const ElementType & type) { return 1; }

  virtual bool checkHomogeneity(const ByElementTypeArray<T> & field_to_check,
                                UInt & nb_comp, UInt & nb_elem) {
    type_iterator tit;
    type_iterator end;

    tit = field.firstType(spatial_dimension, ghost_type, element_kind);
    end = field.lastType(spatial_dimension, ghost_type, element_kind);

    nb_elem = 0;
    nb_comp = 0;
    UInt nb_data = getNbDataPerElem(*tit);

    bool homogen = true;
    if(tit != end) {
      nb_comp = field_to_check(*tit, ghost_type).getNbComponent() * nb_data;
      for(;tit != end; ++tit) {
        const Array<T> & vect = field_to_check(*tit, ghost_type);
        UInt nb_data_cur = getNbDataPerElem(*tit);
        UInt nb_element = filtered ? (*filter)(*tit, ghost_type).getSize() : vect.getSize()/ nb_data_cur;
        UInt nb_comp_cur = vect.getNbComponent() * nb_data_cur;
        if(homogen && nb_comp != nb_comp_cur) homogen = false;
        nb_elem += nb_element;
      }

      if(!homogen) nb_comp = 0;
    }

    return homogen;
  }

protected:
  const ByElementTypeArray<T> & field;
  UInt nb_total_element;
  UInt spatial_dimension;
  GhostType ghost_type;
  ElementKind element_kind;
  bool homogeneous;
  UInt n, itn, out_n;
  const ByElementTypeArray<UInt> * filter;
};

/* -------------------------------------------------------------------------- */
template<typename T, class iterator_type, template<class> class ret_type>
class DumperIOHelper::GenericElementalField<T, iterator_type, ret_type, true> : public Field {
public:
  GenericElementalField(const ByElementTypeArray<T> & field,
                        UInt spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_not_defined,
                        const ByElementTypeArray<UInt> * filter = NULL) :
    field(field), spatial_dimension(spatial_dimension),
    ghost_type(ghost_type), element_kind(element_kind), n(0), itn(0), filter(filter) {
    AKANTU_DEBUG_ASSERT(filter != NULL , "Filter problem!");
    UInt nb_component;
    homogeneous = checkHomogeneity(field, nb_component, nb_total_element);
    if(homogeneous && n == 0) this->n = nb_component;
  }

  GenericElementalField(const ByElementTypeArray<T> & field,
                        UInt n,
                        UInt spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_not_defined,
                        const ByElementTypeArray<UInt> * filter = NULL) :
    field(field), spatial_dimension(spatial_dimension),
    ghost_type(ghost_type), element_kind(element_kind), n(n), itn(0), filter(filter) {
    AKANTU_DEBUG_ASSERT(filter != NULL , "Filter problem!");
    UInt nb_component;
    homogeneous = checkHomogeneity(field, nb_component, nb_total_element);
  }

  typedef iterator_type iterator;
  typedef typename iterator::internal_iterator internal_iterator;
  typedef typename ByElementTypeArray<UInt>::type_iterator type_iterator;
  /* ------------------------------------------------------------------------ */
  virtual iterator begin() {
    type_iterator tit;
    type_iterator end;

    tit = filter->firstType(spatial_dimension, ghost_type, element_kind);
    end = filter->lastType(spatial_dimension, ghost_type, element_kind);

    const Array<T> & vect = field(*tit, ghost_type);
    const Array<UInt> & typed_filter = (*filter)(*tit, ghost_type);
    UInt * fit = typed_filter.storage();

    UInt nb_data = getNbDataPerElem(*tit);
    UInt nb_component = vect.getNbComponent();
    UInt size = vect.getSize() / nb_data;
    UInt ln = nb_component;
    if(itn != 0) ln = itn;

    internal_iterator it = iterator_helper<T, ret_type>::begin(vect, ln, nb_component / ln * nb_data, size);
    iterator rit = iterator(field, n, tit, end, it, *tit, ghost_type, filter, fit);
    rit.setPadding(padding_n, padding_m);
    return rit;
  }

  virtual iterator end  () {
    type_iterator tit;
    type_iterator end;

    tit = filter->firstType(spatial_dimension, ghost_type, element_kind);
    end = filter->lastType(spatial_dimension, ghost_type, element_kind);

    ElementType type = *tit;
    for (; tit != end; ++tit) type = *tit;

    const Array<T> & vect = field(type, ghost_type);
    UInt * fit = NULL;
    const Array<UInt> & typed_filter = (*filter)(type, ghost_type);
    fit = typed_filter.storage()+typed_filter.getSize();

    UInt nb_data = getNbDataPerElem(type);
    UInt nb_component = vect.getNbComponent();
    UInt size = vect.getSize() / nb_data;
    UInt ln = nb_component;
    if(itn != 0 ) ln = itn;
    internal_iterator it = iterator_helper<T, ret_type>::end(vect, ln, nb_component / ln * nb_data, size);
    iterator rit =  iterator(field, n, end, end, it, type, ghost_type, filter, fit);
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
  UInt size() {
    UInt nb_component;
    checkHomogeneity(field, nb_component, nb_total_element);
    return nb_total_element;
  }

  iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }

protected:
  virtual UInt getNbDataPerElem(__attribute__((unused)) const ElementType & type) { return 1; }

  virtual bool checkHomogeneity(const ByElementTypeArray<T> & field_to_check,
                                UInt & nb_comp, UInt & nb_elem) {
    type_iterator tit;
    type_iterator end;

    tit = filter->firstType(spatial_dimension, ghost_type, element_kind);
    end = filter->lastType(spatial_dimension, ghost_type, element_kind);

    nb_elem = 0;
    nb_comp = 0;
    UInt nb_data = getNbDataPerElem(*tit);

    bool homogen = true;
    if(tit != end) {
      nb_comp = field_to_check(*tit, ghost_type).getNbComponent() * nb_data;
      for(;tit != end; ++tit) {
        const Array<T> & vect = field_to_check(*tit, ghost_type);
        UInt nb_data_cur = getNbDataPerElem(*tit);
        UInt nb_element = (*filter)(*tit, ghost_type).getSize();
        UInt nb_comp_cur = vect.getNbComponent() * nb_data_cur;
        if(homogen && nb_comp != nb_comp_cur) homogen = false;
        nb_elem += nb_element;
      }

      if(!homogen) nb_comp = 0;
    }

    return homogen;
  }

protected:
  const ByElementTypeArray<T> & field;
  UInt nb_total_element;
  UInt spatial_dimension;
  GhostType ghost_type;
  ElementKind element_kind;
  bool homogeneous;
  UInt n, itn;
  const ByElementTypeArray<UInt> * filter;
};



/* -------------------------------------------------------------------------- */
template<typename T, template<typename> class ret_type, bool filtered>
class DumperIOHelper::ElementalField : public GenericElementalField<T, elemental_field_iterator<T, ret_type, filtered>, ret_type, filtered> {
public:
  typedef elemental_field_iterator<T, ret_type, filtered> iterator;

  ElementalField(const ByElementTypeArray<T> & field,
                 UInt spatial_dimension = _all_dimensions,
                 GhostType ghost_type = _not_ghost,
                 ElementKind element_kind = _ek_not_defined,
                 const ByElementTypeArray<UInt> * filter = NULL) :
    GenericElementalField<T, iterator, ret_type, filtered>(field, 0, spatial_dimension,
                                                           ghost_type, element_kind, filter) { }
};

template<typename T, bool filtered>
class DumperIOHelper::ElementalField<T, Matrix, filtered> : public GenericElementalField<T, elemental_field_iterator<T, Matrix, filtered>, Matrix, filtered> {
public:
  typedef elemental_field_iterator<T, Matrix, filtered> iterator;

  ElementalField(const ByElementTypeArray<T> & field,
                 UInt n,
                 UInt spatial_dimension = _all_dimensions,
                 GhostType ghost_type = _not_ghost,
                 ElementKind element_kind = _ek_not_defined,
                 const ByElementTypeArray<UInt> * filter = NULL) :
    GenericElementalField<T, iterator, Matrix, filtered>(field, n, spatial_dimension,
                                                         ghost_type, element_kind, filter) { }
};



/* -------------------------------------------------------------------------- */
class DumperIOHelper::FilteredConnectivityField : public GenericElementalField<UInt, filtered_connectivity_field_iterator, Vector, true> {
public:
  typedef filtered_connectivity_field_iterator iterator;
  typedef GenericElementalField<UInt, filtered_connectivity_field_iterator, Vector, true> parent;

  FilteredConnectivityField(const ByElementTypeArray<UInt> & field,
                            const ByElementTypeArray<UInt> & elemental_filter,
                            const Array<UInt> & nodal_filter,
                            UInt spatial_dimension = _all_dimensions,
                            GhostType ghost_type = _not_ghost,
                            ElementKind element_kind = _ek_not_defined) :
    parent(field, 0, spatial_dimension, ghost_type, element_kind, &elemental_filter),
    nodal_filter(nodal_filter) { }

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

private:
  const Array<UInt> & nodal_filter;
};

/* -------------------------------------------------------------------------- */
template<bool filtered>
class DumperIOHelper::ElementTypeField : public GenericElementalField<UInt,
                                                                      element_type_field_iterator<filtered>,
                                                                      Vector, filtered> {
public:
  typedef element_type_field_iterator<filtered> iterator;
private:
  typedef GenericElementalField<UInt, iterator, Vector, filtered> parent;
public:
  /* ------------------------------------------------------------------------ */
  ElementTypeField(const Mesh & mesh,
                   UInt spatial_dimension = _all_dimensions,
                   GhostType ghost_type = _not_ghost,
                   ElementKind element_kind = _ek_not_defined,
                 const ByElementTypeArray<UInt> * filter = NULL) :
    parent(mesh.getConnectivities(), 0, spatial_dimension, ghost_type, element_kind, filter) {
    this->homogeneous = true;
  }

  virtual void registerToDumper(__attribute__((unused)) const std::string & id,
                           iohelper::Dumper & dupmer) {
    dupmer.addElemDataField("element_type", *this);
  }

  UInt getDim() { return 1; }
};


/* -------------------------------------------------------------------------- */
template<bool filtered>
class DumperIOHelper::ElementPartitionField : public GenericElementalField<UInt,
                                                                           element_partition_field_iterator<filtered>,
                                                                           Vector, filtered> {
public:
  typedef element_partition_field_iterator<filtered> iterator;
private:
  typedef GenericElementalField<UInt, iterator, Vector, filtered> parent;
public:
  /* ------------------------------------------------------------------------ */
  ElementPartitionField(const Mesh & mesh,
                        UInt spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_not_defined,
                 const ByElementTypeArray<UInt> * filter = NULL) :
    parent(mesh.getConnectivities(), 0, spatial_dimension, ghost_type, element_kind, filter) {
    this->homogeneous = true;
  }

  UInt getDim() { return 1; }
};
