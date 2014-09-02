/**
 * @file   boundary_inline_impl.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date   Wed Mar 06 09:30:00 2013
 *
 * @brief  Stores information relevent to the notion of domain boundary and surfaces.
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
#include "element_group.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

template <typename T, template <bool> class dump_type>
dumper::Field * GroupManager
::createElementalField(const ElementTypeMapArray<T> & field, 
		       const std::string & group_name,
		       UInt spatial_dimension,
		       const ElementKind & kind,
		       ElementTypeMap<UInt> nb_data_per_elem){
  
  if (&field == NULL) return NULL;
  if (group_name == "all") 
    return this->createElementalField<dump_type<false> >(field,group_name,spatial_dimension,kind,nb_data_per_elem);
  else return this->createElementalFilteredField<dump_type<true> >(field,group_name,spatial_dimension,kind,nb_data_per_elem);
}


/* -------------------------------------------------------------------------- */
template <typename T, template <class> class T2,template <class,template <class> class,bool> class dump_type>
dumper::Field * GroupManager
::createElementalField(const ElementTypeMapArray<T> & field, 
		       const std::string & group_name,
		       UInt spatial_dimension,
		       const ElementKind & kind,
		       ElementTypeMap<UInt> nb_data_per_elem){
  
  if (&field == NULL) return NULL;
  if (group_name == "all") 
    return this->createElementalField<dump_type<T,T2,false> >(field,group_name,spatial_dimension,kind,nb_data_per_elem);
  else 
    return this->createElementalFilteredField<dump_type<T,T2,true>  >(field,group_name,spatial_dimension,kind,nb_data_per_elem);
}


/* -------------------------------------------------------------------------- */

template <typename T, 
	  /// type of InternalMaterialField
	  template<typename T, bool filtered> class dump_type>
dumper::Field * GroupManager::createElementalField(const ElementTypeMapArray<T> & field, 
						   const std::string & group_name,
						   UInt spatial_dimension,
						   const ElementKind & kind,
						   ElementTypeMap<UInt> nb_data_per_elem){

  if (&field == NULL) return NULL;
  if (group_name == "all") 
    return this->createElementalField<dump_type<T,false> >(field,group_name,spatial_dimension,kind,nb_data_per_elem);
  else 
    return this->createElementalFilteredField<dump_type<T,true>  >(field,group_name,spatial_dimension,kind,nb_data_per_elem);
}


/* -------------------------------------------------------------------------- */

template <typename dump_type,typename field_type>
dumper::Field * GroupManager::createElementalField(const field_type & field, 
						   const std::string & group_name,
						   UInt spatial_dimension,
						   const ElementKind & kind,
						   ElementTypeMap<UInt> nb_data_per_elem){
  
  if (&field == NULL) return NULL;
  if (group_name != "all") throw;

  dumper::Field * dumper  = new dump_type(field,spatial_dimension,_not_ghost,kind);
  dumper->setNbDataPerElem(nb_data_per_elem);
  return dumper;
}

/* -------------------------------------------------------------------------- */

template <typename dump_type,typename field_type>
dumper::Field * GroupManager::createElementalFilteredField(const field_type & field, 
							   const std::string & group_name,
							   UInt spatial_dimension,
							   const ElementKind & kind,
							   ElementTypeMap<UInt> nb_data_per_elem){
  
  if (&field == NULL) return NULL;
  if (group_name == "all") throw;

  typedef typename field_type::type T;
  ElementGroup & group = this->getElementGroup(group_name);  
  UInt dim = group.getDimension();
  if (dim != spatial_dimension) throw;
  const ElementTypeMapArray<UInt> & elemental_filter = group.getElements();

  ElementTypeMapArrayFilter<T> * filtered = 
    new ElementTypeMapArrayFilter<T>(field,elemental_filter,nb_data_per_elem);
  
  dumper::Field * dumper  = new dump_type(*filtered,dim,_not_ghost,kind);
  dumper->setNbDataPerElem(nb_data_per_elem);

  return dumper;
}

/* -------------------------------------------------------------------------- */

template <typename type, bool flag, template<class,bool> class ftype>
dumper::Field * GroupManager::createNodalField(const ftype<type,flag> * field,
					       const std::string & group_name,
					       UInt padding_size){

  if (field == NULL) return NULL;
  if (group_name == "all"){
    typedef typename dumper::NodalField<type, false> DumpType;
    DumpType * dumper = new DumpType(*field, 0, 0, NULL);
    dumper->setPadding(padding_size);
    return dumper;
  }
  else {
    ElementGroup & group = this->getElementGroup(group_name);  
    const Array<UInt> * nodal_filter = &(group.getNodes());
    typedef typename dumper::NodalField<type, true> DumpType;
    DumpType * dumper = new DumpType(*field, 0, 0, nodal_filter);
    dumper->setPadding(padding_size);
    return dumper;
  }
  return NULL;
}

/* -------------------------------------------------------------------------- */

template <typename type, bool flag, template<class,bool> class ftype>
dumper::Field * GroupManager::createStridedNodalField(const ftype<type,flag> * field,
						      const std::string & group_name,
						      UInt size, UInt stride, 
						      UInt padding_size){

  if (field == NULL) return NULL;
  if (group_name == "all"){
    typedef typename dumper::NodalField<type, false> DumpType;
    DumpType * dumper = new DumpType(*field, size, stride, NULL);
    dumper->setPadding(padding_size);
    return dumper;
  }
  else {
    ElementGroup & group = this->getElementGroup(group_name);  
    const Array<UInt> * nodal_filter = &(group.getNodes());
    typedef typename dumper::NodalField<type, true> DumpType;
    DumpType * dumper = new DumpType(*field, size, stride, nodal_filter);
    dumper->setPadding(padding_size);
    return dumper;
  }
  return NULL;
}

/* -------------------------------------------------------------------------- */


__END_AKANTU__

