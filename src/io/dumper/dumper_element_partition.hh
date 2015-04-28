/**
 * @file   dumper_element_partition.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  ElementPartition field
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

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__

/* -------------------------------------------------------------------------- */
template<class types>
class element_partition_field_iterator
  : public element_iterator<types, element_partition_field_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef element_iterator<types, dumper::element_partition_field_iterator> parent;
  typedef typename types::return_type return_type;
  typedef typename types::array_iterator array_iterator;
  typedef typename types::field_type field_type;



  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  element_partition_field_iterator(const field_type & field,
				   const typename field_type::type_iterator & t_it,
				   const typename field_type::type_iterator & t_it_end,
				   const array_iterator & array_it,
				   const array_iterator & array_it_end,
				   const GhostType ghost_type = _not_ghost) :
    parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) {
    prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  return_type operator*() {
    return return_type(1, prank);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  UInt prank;
};


/* -------------------------------------------------------------------------- */
template<bool filtered = false>
class ElementPartitionField :
  public GenericElementalField<SingleType<UInt,Vector,filtered>,
			       element_partition_field_iterator> {
public:

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  typedef SingleType<UInt,Vector,filtered> types;
  typedef element_partition_field_iterator<types> iterator;
  typedef GenericElementalField<types,element_partition_field_iterator> parent;
  typedef typename types::field_type field_type;

public:

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  ElementPartitionField(const field_type & field,
			UInt spatial_dimension = _all_dimensions,
			GhostType ghost_type = _not_ghost,
			ElementKind element_kind = _ek_not_defined) :
    parent(field, spatial_dimension, ghost_type, element_kind) {
    this->homogeneous = true;
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  UInt getDim() { return 1; }
};

/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
__END_AKANTU__
