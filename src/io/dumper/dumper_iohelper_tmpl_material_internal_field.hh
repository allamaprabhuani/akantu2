/**
 * @file   dumper_iohelper_tmpl_material_internal_field.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Oct 25 14:41:17 2012
 *
 * @brief  description of material internal field
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
/* Element material field iterator/container                                  */
/* -------------------------------------------------------------------------- */
template<typename T, template<class> class ret_type>
class DumperIOHelper::internal_material_field_iterator : public quadrature_point_iterator< UInt, T, ret_type,
											   internal_material_field_iterator<T, ret_type> > {
public:
  typedef quadrature_point_iterator<UInt, T, ret_type, internal_material_field_iterator<T, ret_type> > parent;
  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
  typedef typename Vector<T>::template const_iterator< ret_type<T> > internal_material_iterator;
public:
  internal_material_field_iterator(const field_type & element_material,
				   UInt n,
				   const typename field_type::type_iterator & t_it,
				   const typename field_type::type_iterator & t_it_end,
				   const internal_iterator & it,
				   ElementType element_type,
				   const GhostType ghost_type = _not_ghost) :
    parent(element_material, 2, t_it, t_it_end,
	   it, element_type, ghost_type), out_n(n) { }

  return_type operator*() {
    const ret_type<UInt> & material_id = *this->vit;
    UInt nb_data = this->fem->getNbQuadraturePoints(*this->tit);
    try {
      const Vector<T> & vect =
	model->getMaterial(material_id(1)).getVector(field_id,
						     *this->tit,
						     this->ghost_type);
      internal_material_iterator it
	= iterator_helper<T, ret_type>::begin(vect, out_n,
					      vect.getNbComponent() / out_n * nb_data,
					      vect.getSize() / nb_data);
      it += material_id(0);

      return *it;
    } catch (...) {
      return return_type();
    }
  }

  void setModel(const SolidMechanicsModel & model) { this->model = &model; }
  void setFieldID(const std::string & field_id) { this->field_id = field_id; }

protected:
  virtual UInt getNbDataPerElem(const ElementType & type) { return 1; }
protected:
  const SolidMechanicsModel * model;
  std::string field_id;
  UInt out_n;
};

/* specialization */
template<typename T>
class DumperIOHelper::internal_material_field_iterator<T, types::Matrix> : public quadrature_point_iterator< UInt, T, types::Matrix,
											   internal_material_field_iterator<T, types::Matrix> > {
public:
  typedef quadrature_point_iterator<UInt, T, types::Matrix, internal_material_field_iterator<T, types::Matrix> > parent;
  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
  typedef typename Vector<T>::template const_iterator< types::Matrix<T> > internal_material_iterator;
public:
  internal_material_field_iterator(const field_type & element_material,
				   UInt n,
				   const typename field_type::type_iterator & t_it,
				   const typename field_type::type_iterator & t_it_end,
				   const internal_iterator & it,
				   ElementType element_type,
				   const GhostType ghost_type = _not_ghost) :
    parent(element_material, 2, t_it, t_it_end,
	   it, element_type, ghost_type), out_n(n) { }

  types::Matrix<T> operator*() {
    const types::Matrix<UInt> & material_id = *this->vit;
    UInt nb_data = this->fem->getNbQuadraturePoints(*this->tit);
    try {
      const Vector<T> & vect =
	model->getMaterial(material_id(1, 0)).getVector(field_id,
							*this->tit,
							this->ghost_type);
      UInt ln = out_n;
      if(out_n == 0) ln = vect.getNbComponent();
      internal_material_iterator it
	= iterator_helper<T, types::Matrix>::begin(vect, ln, vect.getNbComponent() / ln * nb_data,
						   vect.getSize() / nb_data);
      it += material_id(0, 0);

      const types::Matrix<T> & ret = *it;
      if(this->padding_n <= ret.rows() || this->padding_m <= ret.cols())
	return ret;
      else {
	types::Matrix<T> tmp(this->padding_n, this->padding_m * nb_data);
	for (UInt d = 0; d < nb_data; ++d)
	  for (UInt i = 0; i < ret.rows(); ++i)
	    for (UInt j = 0; j < ret.cols(); ++j)
	      tmp(i, j + d * this->padding_m) = ret(i, j + d * ret.cols());
	return tmp;
      }
    } catch (...) {
      return types::Matrix<T>();
    }
  }

  void setModel(const SolidMechanicsModel & model) { this->model = &model; }
  void setFieldID(const std::string & field_id) { this->field_id = field_id; }

protected:
  virtual UInt getNbDataPerElem(const ElementType & type) { return 1; }
protected:
  const SolidMechanicsModel * model;
  std::string field_id;
  UInt out_n;
};

/* -------------------------------------------------------------------------- */
template<typename T, template<class> class ret_type>
class DumperIOHelper::InternalMaterialField : public GenericQuadraturePointsField<UInt,
								  internal_material_field_iterator<T, ret_type>,
								  ret_type> {
public:
  typedef internal_material_field_iterator<T, ret_type> iterator;
private:
  typedef GenericQuadraturePointsField<UInt, iterator, ret_type> parent;
public:
  /* ------------------------------------------------------------------------ */
  InternalMaterialField(const SolidMechanicsModel & model,
			const std::string & field_id,
			UInt spatial_dimension = 0,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    parent(model.getFEM(), model.getElementIndexByMaterial(), spatial_dimension, ghost_type, element_kind),
    model(model), field_id(field_id) {
    init();
  }

  InternalMaterialField(const SolidMechanicsModel & model,
			const std::string & field_id,
			UInt n,
			UInt spatial_dimension = 0,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    parent(model.getFEM(), model.getElementIndexByMaterial(), n, spatial_dimension, ghost_type, element_kind),
    model(model), field_id(field_id) {
    init();
  }

  virtual iterator begin() {
    iterator it = parent::begin();
    it.setItSize(2,1);
    it.setModel(model);
    it.setFieldID(field_id);
    return it;
  }

  virtual iterator end  () {
    iterator it = parent::end();
    it.setItSize(2,1);
    it.setModel(model);
    it.setFieldID(field_id);
    return it;
  }

protected:
  void init() {
    typename ByElementTypeVector<UInt>::type_iterator tit =
      this->field.firstType(this->spatial_dimension, this->ghost_type, this->element_kind);
    typename ByElementTypeVector<UInt>::type_iterator end =
      this->field.lastType (this->spatial_dimension, this->ghost_type, this->element_kind);

    UInt nb_materials = model.getNbMaterials();
    bool found = false;
    for(;tit != end; ++tit) {
      //      UInt nb_quad = this->fem.getNbQuadraturePoints(*tit);
      for (UInt ma = 0; ma < nb_materials; ++ma) {
	const Material & mat = model.getMaterial(ma);
	try {
	  const Vector<Real> & vect __attribute__((unused)) = mat.getVector(field_id, *tit, this->ghost_type);
	  found = true;
	} catch (...) { };
      }
    }

    if(!found) AKANTU_DEBUG_ERROR("The field " << field_id << " does not exists in any materials");
  }

  virtual UInt getNbDataPerElem(const ElementType & type) { return 1; }

protected:
  const SolidMechanicsModel & model;
  std::string field_id;
};
