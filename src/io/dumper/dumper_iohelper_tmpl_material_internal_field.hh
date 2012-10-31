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
/* -------------------------------------------------------------------------- */
template<class T, template<class> class R>
class DumperIOHelper::MaterialPaddingHelper : public PaddingHelper<T, R> {
public:
  MaterialPaddingHelper(const SolidMechanicsModel & model) : model(model) { }

  inline R<T> pad(const R<T> & in, UInt padding_m, UInt padding_n, UInt nb_data, UInt material_id) {
    return PaddingHelper<T, R>::pad(in, padding_m, padding_n, nb_data);
  }

protected:
  const SolidMechanicsModel & model;
};

/* -------------------------------------------------------------------------- */
template<class T, template<class> class R>
class DumperIOHelper::StressPaddingHelper : public MaterialPaddingHelper<T, R> {
public:
  StressPaddingHelper(const SolidMechanicsModel & model) : MaterialPaddingHelper<T, R>(model) {
    AKANTU_DEBUG_ERROR("This class exists only to pad stress (types::Matrix<Real>) to 3D");
  }

  inline R<T> pad(const R<T> & in, UInt padding_m, UInt padding_n, UInt nb_data, UInt material_id) {}
protected:
  const Material & material;
};


template<>
class DumperIOHelper::StressPaddingHelper<Real,
					  types::Matrix> : public MaterialPaddingHelper<Real,
											types::Matrix> {
public:
  StressPaddingHelper(const SolidMechanicsModel & model) :
  MaterialPaddingHelper<Real, types::Matrix>(model) { }

  inline types::Matrix<Real> pad(const types::Matrix<Real> & in,
				 UInt padding_m, UInt padding_n, UInt nb_data,
				 UInt material_id) {
    if(padding_m <= in.rows() && padding_n * nb_data <= in.cols())
      return in;
    else {
      AKANTU_DEBUG_ASSERT(padding_m == 3 && padding_n == 3, "This function can only pad to 3D");
      types::Matrix<Real> stress =
	MaterialPaddingHelper<Real, types::Matrix>::pad(in, 3, 3, nb_data, material_id);
      if(in.rows() == 2 && in.cols() == 2 * nb_data) {
	const Material & material = model.getMaterial(material_id);
	bool plane_strain = !material.getParam<bool>("Plane_Stress");
	if(plane_strain) {
	  Real nu = material.getParam<Real>("nu");
	  for (UInt d = 0; d < nb_data; ++d) {
	    stress(2, 2 + 3*d) = nu * (stress(0, 0 + 3*d) + stress(1, 1 + 3*d));
	  }
	}
      }
      return stress;
    }
  }
};

/* -------------------------------------------------------------------------- */
template<class T, template<class> class R>
class DumperIOHelper::StrainPaddingHelper : public MaterialPaddingHelper<T, R> {
public:
  StrainPaddingHelper(const SolidMechanicsModel & model) : MaterialPaddingHelper<T, R>(model) {
    AKANTU_DEBUG_ERROR("This class exists only to pad strain (types::Matrix<Real>) to 3D");
  }

  inline R<T> pad(const R<T> & in, UInt padding_m, UInt padding_n, UInt nb_data, UInt material_id) {}
};


template<>
class DumperIOHelper::StrainPaddingHelper<Real, types::Matrix> : public MaterialPaddingHelper<Real, types::Matrix> {
public:
  StrainPaddingHelper(const SolidMechanicsModel & model) :
  MaterialPaddingHelper<Real, types::Matrix>(model) { }

  inline types::Matrix<Real> pad(const types::Matrix<Real> & in,
				 UInt padding_m, UInt padding_n, UInt nb_data,
				 UInt material_id) {
    if(padding_m <= in.rows() && padding_n * nb_data <= in.cols())
      return in;
    else {
      AKANTU_DEBUG_ASSERT(padding_m == 3 && padding_n == 3, "This function can only pad to 3D");
      types::Matrix<Real> strain =
	MaterialPaddingHelper<Real, types::Matrix>::pad(in, 3, 3, nb_data, material_id);
      if(in.rows() == 2 && in.cols() == 2 * nb_data) {
	const Material & material = model.getMaterial(material_id);
	bool plane_stress = material.getParam<bool>("Plane_Stress");
	if(plane_stress) {
	  Real nu = material.getParam<Real>("nu");
	  for (UInt d = 0; d < nb_data; ++d) {
	    strain(2, 2 + 3*d) = nu / (nu - 1) * (strain(0, 0 + 3*d) + strain(1, 1 + 3*d));
	  }
	}
      }
      return strain;
    }
  }
};


/* -------------------------------------------------------------------------- */
/* Element material field iterator/container                                  */
/* -------------------------------------------------------------------------- */
template<typename T,
	 template<class> class ret_type,
	 template<typename, template<class> class> class padding_helper_type,
	 template<typename, template<class> class> class int_iterator>
class DumperIOHelper::generic_internal_material_field_iterator : public quadrature_point_iterator< UInt, T, ret_type,
												   int_iterator<T, ret_type> > {
public:
  typedef quadrature_point_iterator<UInt, T, ret_type, int_iterator<T, ret_type> > parent;
  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
  typedef typename Vector<T>::template const_iterator< ret_type<T> > internal_material_iterator;
public:
  generic_internal_material_field_iterator(const field_type & element_material,
					   UInt n,
					   const typename field_type::type_iterator & t_it,
					   const typename field_type::type_iterator & t_it_end,
					   const internal_iterator & it,
					   ElementType element_type,
					   const GhostType ghost_type = _not_ghost) :
    parent(element_material, 2, t_it, t_it_end,
	   it, element_type, ghost_type),
    out_n(n), model(NULL), padding_helper(NULL) { }

  ~generic_internal_material_field_iterator() { delete padding_helper; }

  generic_internal_material_field_iterator(const generic_internal_material_field_iterator & it) : parent(it) {
    out_n = it.out_n;
    field_id = it.field_id;
    if(it.model) {
      model = it.model;
      padding_helper = new padding_helper_type<T, ret_type>(*model);
    }
  }

  return_type operator*() {
    const ret_type<UInt> & material_id = *this->vit;
    UInt nb_data = this->fem->getNbQuadraturePoints(*this->tit);
    try {
      const Vector<T> & vect =
	model->getMaterial(material_id[1]).getVector(field_id,
						     *this->tit,
						     this->ghost_type);
      internal_material_iterator it
	= iterator_helper<T, ret_type>::begin(vect, out_n,
					      vect.getNbComponent() / out_n * nb_data,
					      vect.getSize() / nb_data);
      it += material_id[0];

      return padding_helper->pad(*it, this->padding_m, this->padding_n, nb_data, material_id[1]);
    } catch (...) {
      return return_type();
    }
  }

  void setModel(const SolidMechanicsModel & model) {
    this->model = &model;
    this->padding_helper = new padding_helper_type<T, ret_type>(model);
  }
  void setFieldID(const std::string & field_id) { this->field_id = field_id; }

protected:
  virtual UInt getNbDataPerElem(const ElementType & type) { return 1; }
protected:
  UInt out_n;
  const SolidMechanicsModel * model;
  std::string field_id;
  padding_helper_type<T, ret_type> * padding_helper;
};

/* -------------------------------------------------------------------------- */
template<typename T, template<class> class ret_type>
class DumperIOHelper::internal_material_field_iterator :
  public generic_internal_material_field_iterator<T, ret_type,
						  MaterialPaddingHelper,
						  internal_material_field_iterator> {
public:
  typedef generic_internal_material_field_iterator<T, ret_type, MaterialPaddingHelper,
						   internal_material_field_iterator> parent;
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
    parent(element_material, n, t_it, t_it_end,
	   it, element_type, ghost_type) {  }
};

/* -------------------------------------------------------------------------- */
template<typename T, template<class> class ret_type>
class DumperIOHelper::material_stress_field_iterator :
  public generic_internal_material_field_iterator<T, ret_type, StressPaddingHelper,
						  material_stress_field_iterator> {
public:
  typedef generic_internal_material_field_iterator<T, ret_type, StressPaddingHelper,
						   material_stress_field_iterator> parent;
  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
  typedef typename Vector<T>::template const_iterator< ret_type<T> > internal_material_iterator;
public:
  material_stress_field_iterator(const field_type & element_material,
				 UInt n,
				 const typename field_type::type_iterator & t_it,
				 const typename field_type::type_iterator & t_it_end,
				 const internal_iterator & it,
				 ElementType element_type,
				 const GhostType ghost_type = _not_ghost) :
    parent(element_material, n, t_it, t_it_end,
	   it, element_type, ghost_type) { }
};

/* -------------------------------------------------------------------------- */
template<typename T, template<class> class ret_type>
class DumperIOHelper::material_strain_field_iterator :
  public generic_internal_material_field_iterator<T, ret_type, StrainPaddingHelper,
						  material_strain_field_iterator> {
public:
  typedef generic_internal_material_field_iterator<T, ret_type, StrainPaddingHelper,
						   material_strain_field_iterator> parent;
  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
  typedef typename Vector<T>::template const_iterator< ret_type<T> > internal_material_iterator;
public:
  material_strain_field_iterator(const field_type & element_material,
				 UInt n,
				 const typename field_type::type_iterator & t_it,
				 const typename field_type::type_iterator & t_it_end,
				 const internal_iterator & it,
				 ElementType element_type,
				 const GhostType ghost_type = _not_ghost) :
    parent(element_material, n, t_it, t_it_end,
	   it, element_type, ghost_type) { }
};

/* -------------------------------------------------------------------------- */
template<typename T, template<class> class ret_type,
	 template<typename, template<class> class> class iterator_type>
class DumperIOHelper::InternalMaterialField : public GenericQuadraturePointsField<UInt,
								  iterator_type<T, ret_type>,
								  ret_type> {
public:
  typedef iterator_type<T, ret_type> iterator;
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
