/**
 * @file   dumper_iohelper_tmpl_homogenizing_field.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Oct 25 14:42:24 2012
 *
 * @brief  description of field homogenizing field
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

#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
/* Homogenizing containers                                                    */
/* -------------------------------------------------------------------------- */
template<typename T, class Container, template<class> class sub_type>
class DumperIOHelper::PaddingHomogenizingFunctor {
public:
  PaddingHomogenizingFunctor(Container & cont) : cont(cont){
  }


  UInt getNbComponent() {
    typename Container::iterator it  = cont.begin();
    typename Container::iterator end = cont.end();
    UInt nb_comp = 0;
    for (; it != end; ++it) nb_comp = std::max(nb_comp, (*it).size());
    return nb_comp;
  }

  sub_type<T> operator()(const sub_type<T> & vect, __attribute__((unused)) const ElementType & type) {
    return vect;
  }

protected:
  Container & cont;
};

/* -------------------------------------------------------------------------- */
template<typename T, class Container, template<class> class sub_type>
class DumperIOHelper::AvgHomogenizingFunctor {
public:
  AvgHomogenizingFunctor(Container & cont) : cont(cont){
  }

  UInt getNbComponent() {
    typename Container::iterator it  = cont.begin();
    typename Container::iterator end = cont.end();
    UInt nb_comp = 0;
    UInt i = 0;
    for (; it != end; ++it) {
      ElementType type = it.getType();
      UInt nb_quad = cont.getFEM().getNbQuadraturePoints(type);
      nb_comp = std::max(nb_comp, (*it).size() / nb_quad);
      ++i;
    }
    return nb_comp;
  }

  sub_type<T> operator()(const sub_type<T> & vect, const ElementType & type) {
    UInt nb_quad = cont.getFEM().getNbQuadraturePoints(type);
    if(nb_quad == 1) {
      return vect;
    } else {
      UInt s = vect.size() / nb_quad;
      sub_type<T> v(s);
      for (UInt i = 0; i < s; ++i) {
	for (UInt q = 0; q < nb_quad; ++q) {
	  v[i] += vect[i + q*s] / Real(nb_quad);
	}
      }
      return v;
    }
  }

protected:
  Container & cont;
};

/* specialization */
template<typename T, class Container>
class DumperIOHelper::AvgHomogenizingFunctor<T, Container, types::Matrix> {
public:
  AvgHomogenizingFunctor(Container & cont) : cont(cont){
  }

  UInt getNbComponent() {
    typename Container::iterator it  = cont.begin();
    typename Container::iterator end = cont.end();
    UInt nb_comp = 0;
    UInt i = 0;
    for (; it != end; ++it) {
      ElementType type = it.getType();
      UInt nb_quad = cont.getFEM().getNbQuadraturePoints(type);
      nb_comp = std::max(nb_comp, (*it).size() / nb_quad);
      ++i;
    }
    return nb_comp;
  }

  types::Matrix<T> operator()(const types::Matrix<T> & vect,
			      const ElementType & type) {
    UInt nb_quad = cont.getFEM().getNbQuadraturePoints(type);
    if(nb_quad == 1) {
      return vect;
    } else {
      UInt n = vect.rows();
      UInt m = vect.cols() / nb_quad;
      types::Matrix<T> v(n, m);
      for (UInt i = 0; i < n; ++i)
	for (UInt j = 0; j < m; ++j)
	  for (UInt q = 0; q < nb_quad; ++q)
	    v(i, j) += vect(i, j + q*m) / Real(nb_quad);
      return v;
    }
  }

protected:
  Container & cont;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template<typename T, template< typename, template<class> class> class Container,
	 template<typename, class, template<class> class> class Funct,
	 template<typename> class ret_type>
class DumperIOHelper::HomogenizedField : public Field {
protected:
  typedef typename Container<T, ret_type>::iterator sub_iterator;
  typedef Funct<T, Container<T, ret_type>, ret_type> Functor;
public:
  class iterator {
  public:
    iterator(const sub_iterator & it, Functor & func) : it(it), func(func) {}

    bool operator!=(const iterator & it) { return it.it != this->it; }
    iterator operator++() { ++this->it; return *this; }

    ret_type<T> operator*() {
      return func(*it, it.getType());
    }

  protected:
    sub_iterator it;
    Functor & func;
  };

  HomogenizedField(const SolidMechanicsModel & model,
		   const std::string & field_id,
		   UInt spatial_dimension = 0,
		   GhostType ghost_type = _not_ghost,
		   ElementKind element_kind = _ek_not_defined) :
    cont(model, field_id, spatial_dimension, ghost_type, element_kind),
    funct(cont) {
    nb_component = funct.getNbComponent();
  }

  HomogenizedField(const SolidMechanicsModel & model,
		   const std::string & field_id,
		   UInt n,
		   UInt spatial_dimension = 0,
		   GhostType ghost_type = _not_ghost,
		   ElementKind element_kind = _ek_not_defined) :
    cont(model, field_id, n, spatial_dimension, ghost_type, element_kind),
    funct(cont) {
    nb_component = funct.getNbComponent();
  }

  iterator begin() { return iterator(cont.begin(), funct); }
  iterator end  () { return iterator(cont.end(),   funct); }

  virtual void registerToDumper(const std::string & id,
				iohelper::Dumper & dupmer) {
    dupmer.addElemDataField(id, *this);
  }

  UInt getDim() { 
    if(padding_n && padding_m)
      return padding_m*padding_n;
    else return nb_component;
  }

  UInt size() { return cont.size(); }

  UInt isHomogeneous() { return true; }

  void setPadding(UInt n, UInt m = 1) {
    Field::setPadding(n, m);
    cont.setPadding(n, m);
  }
protected:
  Container<T, ret_type> cont;
  Functor funct;
  UInt nb_component;
};
