/**
 * @file   by_element_type.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 31 11:09:48 2011
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


// Forward declaration of the ByElementTypeVector class
template<class Stored, typename SupportType = ElementType>
class ByElementType;

/* -------------------------------------------------------------------------- */
/* ByElementTypeBase                                                          */
/* -------------------------------------------------------------------------- */
/// Common non templated base class for the ByElementType class
class ByElementTypeBase {
public:
  virtual ~ByElementTypeBase() {};
};

/* -------------------------------------------------------------------------- */
/* ByElementType                                                              */
/* -------------------------------------------------------------------------- */

template<class Stored, typename SupportType>
class ByElementType : public ByElementTypeBase {
private:
  void operator=(const ByElementType &) {};
public:
  ByElementType(const ID & id = "by_element_type",
		const ID & parent_id = "");
  ~ByElementType();

  inline static std::string printType(const SupportType & type, const GhostType & ghost_type);

  inline bool exists(const SupportType & type, const GhostType & ghost_type = _not_ghost) const;

  inline const Stored & operator()(const SupportType & type,
				   const GhostType & ghost_type = _not_ghost) const;

  inline Stored & operator()(const SupportType & type,
			     const GhostType & ghost_type = _not_ghost);

  inline Stored & operator()(const Stored & insert,
			     const SupportType & type,
			     const GhostType & ghost_type = _not_ghost);

  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  typedef std::map<SupportType, Stored> DataMap;
  class type_iterator : private std::iterator<std::forward_iterator_tag, const SupportType> {
  public:
    typedef const SupportType   value_type;
    typedef const SupportType*  pointer;
    typedef const SupportType&  reference;
  protected:
    typedef typename ByElementType<Stored>::DataMap::const_iterator DataMapIterator;
  public:
    type_iterator(DataMapIterator & list_begin,
		  DataMapIterator & list_end,
		  UInt dim,
		  ElementKind ek);

    type_iterator(const type_iterator & it);
    type_iterator() {}

    inline reference operator*();
    inline reference operator*() const;
    inline type_iterator & operator++();
    type_iterator operator++(int);
    inline bool operator==(const type_iterator & other) const;
    inline bool operator!=(const type_iterator & other) const;
    type_iterator & operator=(const type_iterator & other);

  private:
    DataMapIterator list_begin;
    DataMapIterator list_end;
    UInt dim;
    ElementKind kind;
  };

  inline type_iterator firstType(UInt dim = _all_dimensions,
				 GhostType ghost_type = _not_ghost,
				 ElementKind kind = _ek_not_defined) const;
  inline type_iterator lastType(UInt dim = _all_dimensions,
				GhostType ghost_type = _not_ghost,
				ElementKind kind = _ek_not_defined) const;


protected:
  inline DataMap & getData(GhostType ghost_type);
  inline const DataMap & getData(GhostType ghost_type) const;

public:
  AKANTU_GET_MACRO(ID, id, ID);
  AKANTU_SET_MACRO(ID, id, ID);
/* -------------------------------------------------------------------------- */
protected:
  ID id;

  DataMap data;
  DataMap ghost_data;
};

/* -------------------------------------------------------------------------- */
/* Some typedefs                                                              */
/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType = ElementType>
class ByElementTypeArray : public ByElementType<Array<T> *, SupportType>, protected Memory {
protected:
  typedef ByElementType<Array<T> *, SupportType> parent;
  typedef typename parent::DataMap DataMap;
private:
private:
  void operator=(const ByElementType<T, SupportType> &) {};

public:
  typedef typename parent::type_iterator type_iterator;

  ByElementTypeArray() {};

  ByElementTypeArray(const ID & id, const ID & parent_id,
		     const MemoryID & memory_id = 0) :
    parent(id, parent_id), Memory(memory_id) {};

  inline Array<T> & alloc(UInt size,
			   UInt nb_component,
			   const SupportType & type,
			   const GhostType & ghost_type);

  inline void alloc(UInt size,
		    UInt nb_component,
		    const SupportType & type);

  inline const Array<T> & operator()(const SupportType & type,
				      const GhostType & ghost_type = _not_ghost) const;

  inline Array<T> & operator()(const SupportType & type,
				const GhostType & ghost_type = _not_ghost);

  inline void setArray(const SupportType & type,
			const GhostType & ghost_type,
			const Array<T> & vect);

  inline void free();

  inline void onElementsRemoved(const ByElementTypeArray<UInt> & new_numbering);

  virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  ByElementTypeArray operator=(__attribute__((unused)) const ByElementTypeArray & other) {};
};

/// to store data Array<Real> by element type
typedef ByElementTypeArray<Real> ByElementTypeReal;
/// to store data Array<Int> by element type
typedef ByElementTypeArray<Int>  ByElementTypeInt;
/// to store data Array<UInt> by element type
typedef ByElementTypeArray<UInt, ElementType> ByElementTypeUInt;

/// Map of data of type UInt stored in a mesh
typedef std::map<std::string, Array<UInt> *> UIntDataMap;
typedef ByElementType<UIntDataMap, ElementType> ByElementTypeUIntDataMap;

__END_AKANTU__

#endif /* __AKANTU_BY_ELEMENT_TYPE_HH__ */
