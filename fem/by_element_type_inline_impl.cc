/**
 * @file   by_element_type_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug  4 14:41:29 2011
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

/* -------------------------------------------------------------------------- */
/* ByElementType                                                              */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<class Stored>
inline std::string ByElementType<Stored>::printType(const ElementType & type,
						    const GhostType & ghost_type) {
  std::stringstream sstr; sstr << "(" << ghost_type << ":" << type << ")";
  return sstr.str();
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline bool ByElementType<Stored>::exists(ElementType type, GhostType ghost_type) const {
  return this->getData(ghost_type).find(type) != this->getData(ghost_type).end();
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline const Stored & ByElementType<Stored>::operator()(const ElementType & type,
							const GhostType & ghost_type) const {
  typename ByElementType<Stored>::DataMap::const_iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end())
    AKANTU_EXCEPTION("No element of type "
		     << ByElementType<Stored>::printType(type, ghost_type)
		     << " in this ByElementType<"
		     << debug::demangle(typeid(Stored).name()) << "> class");
  return it->second;
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline Stored & ByElementType<Stored>::operator()(const ElementType & type,
						  const GhostType & ghost_type) {
  typename ByElementType<Stored>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  // if(it == this->getData(ghost_type).end())
  //   AKANTU_EXCEPTION("No element of type "
  // 		     << ByElementType<Stored>::printType(type, ghost_type)
  // 		     << " in this ByElementType<"
  // 		     << debug::demangle(typeid(Stored).name()) << "> class");

  if(it == this->getData(ghost_type).end()) {
    ByElementType<Stored>::DataMap & data = this->getData(ghost_type);
    const std::pair<typename DataMap::iterator, bool> & res =
      data.insert(std::pair<ElementType, Stored>(type, Stored()));
    it = res.first;
  }
  return it->second;
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline Stored & ByElementType<Stored>::operator()(const Stored & insert,
						  const ElementType & type,
						  const GhostType & ghost_type) {
  typename ByElementType<Stored>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(it != this->getData(ghost_type).end()) {
    AKANTU_EXCEPTION("Element of type "
  		     << ByElementType<Stored>::printType(type, ghost_type)
   		     << " already in this ByElementType<"
   		     << debug::demangle(typeid(Stored).name()) << "> class");
  } else {
    ByElementType<Stored>::DataMap & data = this->getData(ghost_type);
    const std::pair<typename DataMap::iterator, bool> & res =
      data.insert(std::pair<ElementType, Stored>(type, insert));
    it = res.first;
  }

  return it->second;
}


/* -------------------------------------------------------------------------- */
template<class Stored>
inline typename ByElementType<Stored>::DataMap & ByElementType<Stored>::getData(GhostType ghost_type) {
  if(ghost_type == _not_ghost) return data;
  else return ghost_data;
}

/* -------------------------------------------------------------------------- */
template<class Stored>
inline const typename ByElementType<Stored>::DataMap & ByElementType<Stored>::getData(GhostType ghost_type) const {
  if(ghost_type == _not_ghost) return data;
  else return ghost_data;
}

/* -------------------------------------------------------------------------- */
/// Works only if stored is a pointer to a class with a printself method
template<class Stored>
void ByElementType<Stored>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "ByElementType<" << debug::demangle(typeid(Stored).name()) << "> [" << std::endl;
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const ByElementType<Stored>::DataMap & data = getData(gt);
    typename ByElementType<Stored>::DataMap::const_iterator it;
    for(it = data.begin(); it != data.end(); ++it) {
      stream << space << space << ByElementType<Stored>::printType(it->first, gt) << " [" << std::endl;
      it->second->printself(stream, indent + 3);
    }
  }
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template<class Stored>
ByElementType<Stored>::ByElementType(const ID & id,
				     const ID & parent_id) {
  AKANTU_DEBUG_IN();

  std::stringstream sstr;
  if(parent_id != "") sstr << parent_id << ":";
  sstr << id;

  this->id = sstr.str();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class Stored>
ByElementType<Stored>::~ByElementType() {

}

/* -------------------------------------------------------------------------- */
template <typename T>
inline Vector<T> & ByElementTypeVector<T>::alloc(UInt size,
						 UInt nb_component,
						 const ElementType & type,
						 const GhostType & ghost_type) {
  std::string ghost_id = "";
  if (ghost_type == _ghost) ghost_id = ":ghost";

  Vector<T> * tmp;

  typename ByElementTypeVector<T>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end()) {
    std::stringstream sstr; sstr << this->id << ":" << type << ghost_id;
    tmp = &(Memory::alloc<T>(sstr.str(), size,
			    nb_component, 0));
    std::stringstream sstrg; sstrg << ghost_type;
    tmp->setTag(sstrg.str());
    this->getData(ghost_type)[type] = tmp;
  } else {
    AKANTU_DEBUG_WARNING("The vector " << this->id << this->printType(type, ghost_type)
			 << " already exists, it is resized instead of allocated.");
    tmp = it->second;
    it->second->resize(size);
  }

  return *tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ByElementTypeVector<T>::alloc(UInt size,
					  UInt nb_component,
					  const ElementType & type) {
  this->alloc(size, nb_component, type, _not_ghost);
  this->alloc(size, nb_component, type, _ghost);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ByElementTypeVector<T>::free() {
  AKANTU_DEBUG_IN();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const DataMap & data = this->getData(gt);
    typename DataMap::const_iterator it;
    for(it = data.begin(); it != data.end(); ++it) {
      dealloc(it->second->getID());
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const Vector<T> & ByElementTypeVector<T>::operator()(const ElementType & type,
							    const GhostType & ghost_type) const {
  typename ByElementTypeVector<T>::DataMap::const_iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end())
    AKANTU_EXCEPTION("No element of type "
		     << ByElementTypeVector<T>::printType(type, ghost_type)
		     << " in this ByElementTypeVector<"
		     << debug::demangle(typeid(T).name()) << "> class");

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline Vector<T> & ByElementTypeVector<T>::operator()(const ElementType & type,
						      const GhostType & ghost_type) {
  typename ByElementTypeVector<T>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(it == this->getData(ghost_type).end())
    AKANTU_EXCEPTION("No element of type "
		     << ByElementTypeVector<T>::printType(type, ghost_type)
		     << " in this ByElementTypeVector<"
		     << debug::demangle(typeid(T).name()) << "> class (\""
		     << this->id << "\")");

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ByElementTypeVector<T>::setVector(const ElementType & type,
					      const GhostType & ghost_type,
					      const Vector<T> & vect) {
  typename ByElementTypeVector<T>::DataMap::iterator it =
    this->getData(ghost_type).find(type);

  if(AKANTU_DEBUG_TEST(dblWarning) && it != this->getData(ghost_type).end()) {
    AKANTU_DEBUG_WARNING("The Vector " << this->printType(type, ghost_type)
			 << " is already registred, this call can lead to a memory leek.");
  }

  this->getData(ghost_type)[type] = &(const_cast<Vector<T> &>(vect));
}

