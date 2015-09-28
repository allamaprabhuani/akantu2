/**
 * @file   dumper_internal_material_field.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  description of material internal field
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


#ifndef __AKANTU_DUMPER_INTERNAL_MATERIAL_FIELD_HH__
#define __AKANTU_DUMPER_INTERNAL_MATERIAL_FIELD_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_quadrature_point_iterator.hh"
#ifdef AKANTU_IGFEM
#  include "dumper_igfem_material_internal_field.hh"
#endif
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */

template<typename T, bool filtered = false>
class InternalMaterialField
  : public GenericElementalField<SingleType<T,Vector,filtered>,
                                 quadrature_point_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:

  typedef SingleType<T,Vector,filtered> types;
  typedef GenericElementalField<types,quadrature_point_iterator> parent;
  typedef typename types::field_type field_type;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */


  InternalMaterialField(const field_type & field,
                        UInt spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_not_defined) :
    parent(field, spatial_dimension, ghost_type, element_kind){}


};


__END_AKANTU_DUMPER__
__END_AKANTU__

#endif /* __AKANTU_DUMPER_INTERNAL_MATERIAL_FIELD_HH__ */
