/**
 * @file   dumper_iohelper_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 26 21:52:40 2012
 *
 * @brief  implementation of the DumperIOHelper class
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
class DumperIOHelper::iterator_helper {
  typedef typename Vector<T>::template const_iterator< R<T> > internal_iterator;
public:
  static internal_iterator begin(const Vector<T> & vect, UInt m, UInt n, UInt size) {
    return vect.begin_reinterpret(n*m, size);
  }

  static internal_iterator end(const Vector<T> & vect, UInt m, UInt n, UInt size) {
    return vect.end_reinterpret(n*m, size);
  }
};


template<class T>
class DumperIOHelper::iterator_helper<T, types::Matrix> {
  typedef typename Vector<T>::template const_iterator< types::Matrix<T> > internal_iterator;
public:
  static internal_iterator begin(const Vector<T> & vect, UInt m, UInt n, UInt size) {
    return vect.begin_reinterpret(m, n, size);
  }

  static internal_iterator end(const Vector<T> & vect, UInt m, UInt n, UInt size) {
    return vect.end_reinterpret(m, n, size);
  }
};


/* -------------------------------------------------------------------------- */
template<class T, template<class> class R>
class DumperIOHelper::PaddingHelper {
public:
  static inline R<T> pad(const R<T> & in,
		  __attribute__((unused)) UInt padding_m,
		  __attribute__((unused)) UInt padding_n,
		  __attribute__((unused)) UInt nb_data) {
    return in; // trick due to the fact that IOHelper padds the vectors (avoid a copy of data)
  }
};

template<class T>
class DumperIOHelper::PaddingHelper<T, types::Matrix> {
public:
  static inline types::Matrix<T> pad(const types::Matrix<T> & in, UInt padding_m, UInt padding_n, UInt nb_data) {
    if(padding_m <= in.rows() && padding_n*nb_data <= in.cols())
      return in;
    else {
      types::Matrix<T> ret(padding_m, padding_n * nb_data);
      for (UInt d = 0; d < nb_data; ++d)
	for (UInt i = 0; i < in.rows(); ++i)
	  for (UInt j = 0; j < in.cols() / nb_data; ++j)
	    ret(i, j + d * padding_n) = in(i, j + d * in.cols());
      return ret;
    }
  }
};


#include "dumper_iohelper_tmpl_nodal_field.hh"
#include "dumper_iohelper_tmpl_elemental_field.hh"
#include "dumper_iohelper_tmpl_quadrature_points_field.hh"
