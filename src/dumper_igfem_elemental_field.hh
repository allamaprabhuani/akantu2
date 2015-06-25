/**
 * @file   dumper_igfem_elemental_field.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  description of IGFEM elemental fields
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef __AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH__
#define __AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH__
/* -------------------------------------------------------------------------- */
#include "static_communicator.hh"
#include "dumper_field.hh"
#include "dumper_igfem_generic_elemental_field.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */


template<typename T, template <class> class ret = Vector,bool filtered = false>
class IGFEMElementalField
  : public IGFEMGenericElementalField<SingleType<T,ret,filtered>,
				      igfem_elemental_field_iterator> {

public:

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  typedef SingleType<T,ret,filtered> types;
  typedef typename types::field_type field_type;
  typedef elemental_field_iterator<types> iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  IGFEMElementalField(const field_type & field,
		      UInt spatial_dimension = _all_dimensions,
		      GhostType ghost_type = _not_ghost) :
    IGFEMGenericElementalField<types,igfem_elemental_field_iterator>(field,
								     spatial_dimension,
								     ghost_type) { }
};


__END_AKANTU_DUMPER__
__END_AKANTU__

#endif /* __AKANTU_DUMPER_IGFEM_ELEMENTAL_FIELD_HH__ */
