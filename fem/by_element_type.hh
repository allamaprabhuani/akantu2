/**
 * @file   by_element_type.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug  4 14:40:34 2011
 *
 * @brief  storage class by element type
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

#ifndef __AKANTU_BY_ELEMENT_TYPE_HH__
#define __AKANTU_BY_ELEMENT_TYPE_HH__

#include "aka_common.hh"
#include "aka_vector.hh"
#include "aka_memory.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* ByElementType                                                              */
/* -------------------------------------------------------------------------- */

template<class Stored> class ByElementType {
protected:
  typedef std::map<ElementType, Stored> DataMap;
public:
  ByElementType(const ID & id = "by_element_type",
		const ID & parent_id = "");
  ~ByElementType();

  __aka_inline__ static std::string printType(const ElementType & type, const GhostType & ghost_type);

  __aka_inline__ bool exists(ElementType type, GhostType ghost_type = _not_ghost) const;

  __aka_inline__ const Stored & operator()(const ElementType & type,
				   const GhostType & ghost_type = _not_ghost) const;
  __aka_inline__ Stored & operator()(const ElementType & type,
			     const GhostType & ghost_type = _not_ghost);

  __aka_inline__ Stored & operator()(const Stored & insert,
			     const ElementType & type,
			     const GhostType & ghost_type = _not_ghost);

  void printself(std::ostream & stream, int indent = 0) const;

protected:
  __aka_inline__ DataMap & getData(GhostType ghost_type);
  __aka_inline__ const DataMap & getData(GhostType ghost_type) const;

/* -------------------------------------------------------------------------- */
protected:
  ID id;

  DataMap data;
  DataMap ghost_data;
};


/* -------------------------------------------------------------------------- */
/* Some typedefs                                                              */
/* -------------------------------------------------------------------------- */

template <typename T>
class ByElementTypeVector : public ByElementType<Vector<T> *>, protected Memory {
protected:
  typedef typename ByElementType<Vector<T> *>::DataMap DataMap;
public:
  ByElementTypeVector() {};
  // ByElementTypeVector(const ID & id = "by_element_type_vector",
  // 		      const MemoryID & memory_id = 0) :
  //   ByElementType<Vector<T> *>(id, memory_id) {};
  ByElementTypeVector(const ID & id, const ID & parent_id,
		      const MemoryID & memory_id = 0) :
    ByElementType<Vector<T> *>(id, parent_id), Memory(memory_id) {};

  __aka_inline__ Vector<T> & alloc(UInt size,
			   UInt nb_component,
			   const ElementType & type,
			   const GhostType & ghost_type);

  __aka_inline__ void alloc(UInt size,
		    UInt nb_component,
		    const ElementType & type);

  __aka_inline__ const Vector<T> & operator()(const ElementType & type,
				      const GhostType & ghost_type = _not_ghost) const;

  __aka_inline__ Vector<T> & operator()(const ElementType & type,
				const GhostType & ghost_type = _not_ghost);

  __aka_inline__ void setVector(const ElementType & type,
			const GhostType & ghost_type,
			const Vector<T> & vect);

  __aka_inline__ void free();
};

/// to store data Vector<Real> by element type
typedef ByElementTypeVector<Real> ByElementTypeReal;
/// to store data Vector<Int> by element type
typedef ByElementTypeVector<Int>  ByElementTypeInt;
/// to store data Vector<UInt> by element type
typedef ByElementTypeVector<UInt> ByElementTypeUInt;

/// Map of data of type UInt stored in a mesh
typedef std::map<std::string, Vector<UInt> *> UIntDataMap;
typedef ByElementType<UIntDataMap> ByElementTypeUIntDataMap;

/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "by_element_type_inline_impl.cc"
#endif


__END_AKANTU__

#endif /* __AKANTU_BY_ELEMENT_TYPE_HH__ */
