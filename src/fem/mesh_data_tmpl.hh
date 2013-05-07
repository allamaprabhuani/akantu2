/**
 * @file   mesh_data_tmpl.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date   Thu Feb 27 15:00:00 2013
 *
 * @brief  Stores generic data loaded from the mesh file
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

__END_AKANTU__

#include <iostream>

#define MESH_DATA_GET_TYPE(r, data, type) \
template<> \
inline MeshDataTypeCode MeshData::getTypeCode<BOOST_PP_TUPLE_ELEM(2, 1, type)>() const { \
  return BOOST_PP_TUPLE_ELEM(2, 0, type); \
}

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename T>
inline MeshDataTypeCode MeshData::getTypeCode() const {
  AKANTU_DEBUG_ERROR("Type " << debug::demangle(typeid(T).name()) << "not implemented by MeshData.");
}

/* -------------------------------------------------------------------------- */
BOOST_PP_SEQ_FOR_EACH(MESH_DATA_GET_TYPE, void, AKANTU_MESH_DATA_TYPES)
#undef MESH_DATA_GET_TYPE

/* -------------------------------------------------------------------------- */
inline MeshDataTypeCode MeshData::getTypeCode(const std::string name) const {
  TypeCodeMap::const_iterator it = typecode_map.find(name);
  if(it == typecode_map.end()) AKANTU_EXCEPTION("No dataset named " << name << " found.");
  return it->second;
}

/* -------------------------------------------------------------------------- */
template<typename T>
void MeshData::registerElementalData(const std::string & name) {
  ElementalDataMap::iterator it = elemental_data.find(name);
  if(it == elemental_data.end()) {
    allocElementalData<T>(name);
  } else{
    AKANTU_DEBUG_WARNING("Data named " << name << " already registered.");
  }
}

/* -------------------------------------------------------------------------- */
#define AKANTU_MESH_DATA_CASE_MACRO(r, name, elem)	\
  case BOOST_PP_TUPLE_ELEM(2, 0, elem) : { registerElementalData<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(name); break; }

inline void MeshData::registerElementalData(const std::string & name, MeshDataTypeCode type) {
  switch(type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_MESH_DATA_CASE_MACRO, name, AKANTU_MESH_DATA_TYPES)
  default : AKANTU_DEBUG_ERROR("Type " << type << "not implemented by MeshData.");
  }
}
#undef AKANTU_MESH_DATA_CASE_MACRO
/* -------------------------------------------------------------------------- */
template<typename T>
ByElementTypeArray<T> * MeshData::allocElementalData(const std::string & name) {
  ByElementTypeArray<T> * dataset = new ByElementTypeArray<T>(name, id, memory_id);
  elemental_data[name] = dataset;
  typecode_map[name] = getTypeCode<T>();
  return dataset;
}

/* -------------------------------------------------------------------------- */
template<typename T>
const ByElementTypeArray<T> & MeshData::getElementalData(const std::string & name) const {
  ElementalDataMap::const_iterator it = elemental_data.find(name);
  if(it == elemental_data.end()) AKANTU_EXCEPTION("No dataset named " << name << " found.");
  return dynamic_cast<const ByElementTypeArray<T> &>(*(it->second));
}

/* -------------------------------------------------------------------------- */
template<typename T>
ByElementTypeArray<T> & MeshData::getElementalData(const std::string & name) {
  ElementalDataMap::iterator it = elemental_data.find(name);
  if(it == elemental_data.end()) AKANTU_EXCEPTION("No dataset named " << name << " found.");
  return dynamic_cast<ByElementTypeArray<T> &>(*(it->second));
}

/* -------------------------------------------------------------------------- */
template<typename T>
const Array<T> & MeshData::getElementalDataArray(const std::string & name,
                                                 const ElementType & elem_type,
                                                 const GhostType & ghost_type) const {
  ElementalDataMap::const_iterator it = elemental_data.find(name);
  if(it == elemental_data.end()) {
    AKANTU_EXCEPTION("Data named " << name << " not registered for type: " << elem_type << " - ghost_type:" << ghost_type << "!");
  }
  return dynamic_cast<const ByElementTypeArray<T> &>(*(it->second))(elem_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template<typename T>
Array<T> & MeshData::getElementalDataArray(const std::string & name,
                                           const ElementType & elem_type,
                                           const GhostType & ghost_type) {
  ElementalDataMap::iterator it = elemental_data.find(name);
  if(it == elemental_data.end()) {
    AKANTU_EXCEPTION("Data named " << name << " not registered for type: " << elem_type << " - ghost_type:" << ghost_type << "!");
  }
  return dynamic_cast<ByElementTypeArray<T> &>(*(it->second))(elem_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template<typename T>
Array<T> & MeshData::getElementalDataArrayAlloc(const std::string & name,
                                                  const ElementType & elem_type,
                                                  const GhostType & ghost_type,
                                                  UInt nb_component) {
  ElementalDataMap::iterator it = elemental_data.find(name);
  ByElementTypeArray<T> * dataset;
  if(it == elemental_data.end()) {
    dataset = allocElementalData<T>(name);
  } else {
    dataset = dynamic_cast<ByElementTypeArray<T> *>(it->second);
  }
  AKANTU_DEBUG_ASSERT(getTypeCode<T>() == typecode_map.find(name)->second, "Function getElementalDataArrayAlloc called with the wrong type!");
  if(!(dataset->exists(elem_type, ghost_type))) {
    dataset->alloc(0, nb_component, elem_type, ghost_type);
  }
  return (*dataset)(elem_type, ghost_type);
}


/* -------------------------------------------------------------------------- */
#define AKANTU_MESH_DATA_CASE_MACRO(r, name, elem)	\
  case BOOST_PP_TUPLE_ELEM(2, 0, elem) : { \
    nb_comp = getNbComponentTemplated<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(name, el_type, ghost_type); \
    break; \
  } \

inline UInt MeshData::getNbComponent(const std::string name, const ElementType & el_type, const GhostType & ghost_type) const {
  TypeCodeMap::const_iterator it = typecode_map.find(name);
  UInt nb_comp(0);
  if(it == typecode_map.end()) {
    AKANTU_EXCEPTION("Could not determine the type held in dataset " << name << " for type: " << el_type << " - ghost_type:" << ghost_type << ".");
  }
  MeshDataTypeCode type = it->second;
  switch(type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_MESH_DATA_CASE_MACRO, name, AKANTU_MESH_DATA_TYPES)
  default : AKANTU_DEBUG_ERROR("Could not call the correct instance of getNbComponentTemplated."); break;
  }
  return nb_comp;
}
#undef AKANTU_MESH_DATA_CASE_MACRO

/* -------------------------------------------------------------------------- */
template<typename T>
inline UInt MeshData::getNbComponentTemplated(const std::string name, const ElementType & el_type, const GhostType & ghost_type) const {
  return getElementalDataArray<T>(name, el_type, ghost_type).getNbComponent();
}

/* -------------------------------------------------------------------------- */
#define AKANTU_MESH_DATA_CASE_MACRO(r, name, elem)	\
  case BOOST_PP_TUPLE_ELEM(2, 0, elem) : { \
    ByElementTypeArray<BOOST_PP_TUPLE_ELEM(2, 1, elem)> * dataset; \
    dataset = dynamic_cast< ByElementTypeArray<BOOST_PP_TUPLE_ELEM(2, 1, elem)> * >(it->second); \
    exists = dataset->exists(el_type, ghost_type); \
    break; \
  } \

inline void MeshData::getTagNames(StringVector & tags, const ElementType & el_type, const GhostType & ghost_type) const {
  tags.clear();
  bool exists(false);

  ElementalDataMap::const_iterator it     = elemental_data.begin();
  ElementalDataMap::const_iterator it_end = elemental_data.end();
  for(; it != it_end; ++it) {
    MeshDataTypeCode type = getTypeCode(it->first);
    switch(type) {
      BOOST_PP_SEQ_FOR_EACH(AKANTU_MESH_DATA_CASE_MACRO, , AKANTU_MESH_DATA_TYPES)
      default : AKANTU_DEBUG_ERROR("Could not determine the proper type to (dynamic-)cast."); break;
    }
    if(exists) {
      tags.push_back(it->first);
    }
  }
}
#undef AKANTU_MESH_DATA_CASE_MACRO

/* -------------------------------------------------------------------------- */



