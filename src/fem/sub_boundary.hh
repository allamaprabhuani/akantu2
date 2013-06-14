/**
 * @file   sub_boundary.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
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

#ifndef __AKANTU_SUB_BOUNDARY_HH__
#define __AKANTU_SUB_BOUNDARY_HH__

#include <set>
#include "aka_common.hh"
#include "by_element_type.hh"
#include "dumpable.hh"

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
class SubBoundary : private Memory, public Dumpable {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  // Warning: The NodeList attribute holding the nodes indices has to be sorted!
  //          Other classes depend on it (e.g. BoundaryNodalField in the dumper)
  typedef Array<UInt> NodeList;
  typedef ByElementTypeArray<UInt> ElementList;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SubBoundary(const std::string & name, const std::string & id = "sub_boundary", 
	      const MemoryID & memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Node iterator                                                            */
  /* ------------------------------------------------------------------------ */
  class nodes_const_iterator {
  public:
    inline nodes_const_iterator(const nodes_const_iterator &);

  private:
    inline nodes_const_iterator(const NodeList::const_iterator<UInt> &);

  public:
    inline bool operator==(const nodes_const_iterator &) const;
    inline bool operator!=(const nodes_const_iterator &) const;
    inline nodes_const_iterator & operator=(const nodes_const_iterator &);
    inline nodes_const_iterator & operator++();
    inline nodes_const_iterator operator++(int);
    inline const UInt & operator*() const;

  private:
    friend class SubBoundary;
    NodeList::const_iterator<UInt> iter;
  };

  inline nodes_const_iterator nodes_begin() const;
  inline nodes_const_iterator nodes_end() const;

  /* ------------------------------------------------------------------------ */
  /* Methods and accessors                                                    */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Elements, elements, UInt);
  AKANTU_GET_MACRO(Elements, elements, const ByElementTypeArray<UInt> &);
  AKANTU_GET_MACRO(Nodes, nodes, const Array<UInt> &);

  AKANTU_GET_MACRO(Name, name, std::string);
  AKANTU_GET_MACRO(ID, id, const ID &); 
  inline UInt getNbNodes() const;
  void printself(std::ostream & stream) const;
  inline void registerField(const std::string & dumper_name,
			    const std::string & field_name, 
			    DumperIOHelper::Field * field);

private:
  inline void addNode(UInt node_id);
  inline void addElement(const ElementType & elem_type,
			 UInt elem_id,
			 const GhostType & ghost_type);
  void cleanUpNodeList();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ID id;
  friend class Boundary;
  std::string name;
  NodeList & nodes;
  ElementList elements;

};

#include "sub_boundary_inline_impl.cc"

__END_AKANTU__


#endif /* __AKANTU_SUB_BOUNDARY_HH__ */
