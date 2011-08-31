/**
 * @file   material_non_local_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 25 11:59:39 2011
 *
 * @brief  
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

#include "aka_types.hh"

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
template<typename T>
void MaterialNonLocal::accumulateNeighbours(const ByElementTypeVector<T> & to_accumulate,
					    ByElementTypeVector<T> & accumulated,
					    UInt nb_degree_of_freedom) const {
  std::set< std::pair<ElementType, ElementType> >::const_iterator first_pair_types = existing_pairs.begin();
  std::set< std::pair<ElementType, ElementType> >::const_iterator last_pair_types = existing_pairs.end();

  GhostType ghost_type1, ghost_type2;
  ghost_type1 = ghost_type2 = _not_ghost;

  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Vector<UInt> & pairs =
      pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    const Vector<T> & to_acc = to_accumulate(first_pair_types->second, ghost_type2);
    Vector<T> & acc = accumulated(first_pair_types->first, ghost_type1);

    acc.copy(to_acc);

    Vector<UInt>::const_iterator< types::Vector<UInt> > first_pair = pairs.begin(2);
    Vector<UInt>::const_iterator< types::Vector<UInt> > last_pair  = pairs.end(2);

    typename Vector<T>::template const_iterator< types::Vector<T> > to_acc_it = to_acc.begin(nb_degree_of_freedom);
    typename Vector<T>::template iterator< typename types::Vector<T> > acc_it = acc.begin(nb_degree_of_freedom);

    for(;first_pair != last_pair; ++first_pair) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);
      acc_it[q1] += to_acc_it[q2];
    }
  }
}

