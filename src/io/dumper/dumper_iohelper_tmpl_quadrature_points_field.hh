/**
 * @file   dumper_iohelper_tmpl_quadrature_points_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 26 21:52:40 2012
 *
 * @brief  Description of quadrature points fields
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

#include "fem.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename i_type, typename d_type,
	 template<typename> class ret_type, class daughter>
class DumperIOHelper::generic_quadrature_point_iterator : public element_iterator<i_type, d_type,
									      ret_type, daughter> {
public:
  typedef element_iterator<i_type, d_type, ret_type, daughter> parent;
  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
public:
  generic_quadrature_point_iterator(const field_type & field,
				    UInt n,
				    const typename field_type::type_iterator & t_it,
				    const typename field_type::type_iterator & t_it_end,
				    const internal_iterator & it,
				    ElementType element_type,
				    const GhostType ghost_type = _not_ghost) :
    parent(field, n, t_it, t_it_end, it, element_type, ghost_type) { }

  void setFEM(const FEM & fem) { this->fem = &fem; }

protected:
  virtual UInt getNbDataPerElem(const ElementType & type) { return fem->getNbQuadraturePoints(type); }

protected:
  const FEM * fem;
};

/* -------------------------------------------------------------------------- */
template<typename T, template<typename> class ret_type>
class DumperIOHelper::quadrature_point_iterator : public generic_quadrature_point_iterator<T, T,
											   ret_type, quadrature_point_iterator<T, ret_type> > {
 public:
  typedef generic_quadrature_point_iterator<T, T, ret_type, quadrature_point_iterator> parent;
  typedef typename parent::it_type     it_type;
  typedef typename parent::data_type   data_type;
  typedef typename parent::return_type return_type;
  typedef typename parent::field_type  field_type;
  typedef typename parent::internal_iterator internal_iterator;
public:
  quadrature_point_iterator(const field_type & field,
			    UInt n,
			    const typename field_type::type_iterator & t_it,
			    const typename field_type::type_iterator & t_it_end,
			    const internal_iterator & it,
			    ElementType element_type,
			    const GhostType ghost_type = _not_ghost) :
    parent(field, n, t_it, t_it_end, it, element_type, ghost_type) { }

  return_type operator*() {
    UInt nb_data = this->fem->getNbQuadraturePoints(*this->tit);
    return PaddingHelper<T, ret_type>::pad(*this->vit, this->padding_m, this->padding_n, nb_data);
  }
};

/* -------------------------------------------------------------------------- */
/* Fields type description                                                    */
/* -------------------------------------------------------------------------- */
template<typename T,
	 class iterator_type,
	 template<class> class ret_type>
class DumperIOHelper::GenericQuadraturePointsField : public GenericElementalField<T, iterator_type, ret_type> {
public:
  typedef iterator_type iterator;
  typedef GenericElementalField<T, iterator, ret_type> parent;

  GenericQuadraturePointsField(const FEM & fem,
			       const ByElementTypeVector<T> & field,
			       UInt spatial_dimension = 0,
			       GhostType ghost_type = _not_ghost,
			       ElementKind element_kind = _ek_not_defined) :
    parent(field, spatial_dimension,
	   ghost_type, element_kind), fem(fem) { }

  GenericQuadraturePointsField(const FEM & fem,
			       const ByElementTypeVector<T> & field,
			       UInt n,
			       UInt spatial_dimension = 0,
			       GhostType ghost_type = _not_ghost,
			       ElementKind element_kind = _ek_not_defined) :
    parent(field, n, spatial_dimension,
	   ghost_type, element_kind), fem(fem) { }

  /* ------------------------------------------------------------------------ */
  virtual iterator begin() {
    iterator it = parent::begin();
    it.setFEM(fem);
    return it;
  }

  virtual iterator end  () {
    iterator it = parent::end();
    it.setFEM(fem);
    return it;
  }

  const FEM & getFEM() { return fem; }

protected:
  virtual UInt getNbDataPerElem(const ElementType & type) { return fem.getNbQuadraturePoints(type); }

protected:
  const FEM & fem;
};

/* -------------------------------------------------------------------------- */
template<typename T,
	 template<class> class ret_type,
	 template<typename, template<class> class> class iterator_type>
class DumperIOHelper::QuadraturePointsField : public GenericQuadraturePointsField<T, iterator_type<T, ret_type>, ret_type> {
public:
  typedef iterator_type<T, ret_type> iterator;
  typedef GenericQuadraturePointsField<T, iterator, ret_type> parent;

  QuadraturePointsField(const FEM & fem,
			const ByElementTypeVector<T> & field,
			UInt spatial_dimension = 0,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    parent(fem, field, spatial_dimension,
	   ghost_type, element_kind) { }

  QuadraturePointsField(const FEM & fem,
			const ByElementTypeVector<T> & field,
			UInt n,
			UInt spatial_dimension = 0,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    parent(fem, field, n, spatial_dimension,
	   ghost_type, element_kind) { }
};
