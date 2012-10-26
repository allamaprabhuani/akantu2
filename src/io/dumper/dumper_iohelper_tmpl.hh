/**
 * @file   dumper_iohelper_tmpl.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Oct 24 13:35:18 2012
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
  static internal_iterator begin(const Vector<T> & vect, UInt n, UInt m, UInt size) {
    return vect.begin_reinterpret(n*m, size);
  }

  static internal_iterator end(const Vector<T> & vect, UInt n, UInt m, UInt size) {
    return vect.end_reinterpret(n*m, size);
  }
};


template<class T>
class DumperIOHelper::iterator_helper<T, types::Matrix> {
  typedef typename Vector<T>::template const_iterator< types::Matrix<T> > internal_iterator;
public:
  static internal_iterator begin(const Vector<T> & vect, UInt n, UInt m, UInt size) {
    return vect.begin_reinterpret(n, m, size);
  }

  static internal_iterator end(const Vector<T> & vect, UInt n, UInt m, UInt size) {
    return vect.end_reinterpret(n, m, size);
  }
};


#include "dumper_iohelper_tmpl_nodal_field.hh"

#include "dumper_iohelper_tmpl_elemental_field.hh"
#include "dumper_iohelper_tmpl_quadrature_points_field.hh"

