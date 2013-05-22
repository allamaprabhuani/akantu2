/**
 * @file   boundary.hh
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

#ifndef __AKANTU_BOUNDARY_HH__
#define __AKANTU_BOUNDARY_HH__

#include <set>
#include "aka_common.hh"
#include "by_element_type.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class Mesh;
/* -------------------------------------------------------------------------- */
class SubBoundary;

class Boundary {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::map<std::string, SubBoundary *> BoundaryList;

public:
  typedef std::set<ElementType> BoundaryTypeSet;

  /* ------------------------------------------------------------------------ */
  /* SubBoundary iterator                                                            */
  /* ------------------------------------------------------------------------ */
public:
  template<typename T, typename container_iterator>
  class internal_iterator {

  public:
    inline internal_iterator(const internal_iterator &);

  private:
    explicit inline internal_iterator(const container_iterator &);

  public:
    inline T & operator*() const;
    inline T * operator->() const;
    inline bool operator==(const internal_iterator &) const;
    inline bool operator!=(const internal_iterator &) const;
    inline internal_iterator & operator=(const internal_iterator &);
    inline internal_iterator & operator++();
    inline internal_iterator operator++(int);

  private:
    friend class Boundary;
    container_iterator iter;
  };

  typedef internal_iterator<const SubBoundary,
                            BoundaryList::const_iterator> const_iterator;

  typedef internal_iterator<SubBoundary,
                            BoundaryList::iterator> iterator;

  inline const_iterator begin() const;
  inline const_iterator find(std::string name) const;
  inline const_iterator end() const;

  inline iterator begin();
  inline iterator find(std::string name);
  inline iterator end();

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
   Boundary(const Mesh & mesh, const ID & id = "boundary", const ID & parent_id = "", const MemoryID & memory_id = 0);
   ~Boundary();

  /* ------------------------------------------------------------------------ */
  /* Methods and accessors                                                    */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T> void createBoundariesFromMeshData(const std::string & dataset_name);
  void createBoundariesFromMeshData(const std::string & dataset_name);
  void createBoundariesFromGeometry();

  /// Create a SubBoundary based on a node group
  void createSubBoundaryFromNodeGroup(const std::string & name,
				      const Array<UInt> & node_group);
  inline const SubBoundary & operator()(const std::string & name) const;
  BoundaryTypeSet getBoundaryElementTypes();
  AKANTU_GET_MACRO(NbBoundaries, boundaries.size(), UInt);
  void printself(std::ostream & stream) const;
  void dump();

private:
  void addElementAndNodesToBoundaryAlloc(const std::string & boundary_name,
					 const ElementType & elem_type,
					 UInt elem_id,
					 const GhostType & ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ID id;
  BoundaryList boundaries;
  MemoryID memory_id;
  const Mesh & mesh;
};

#include "boundary_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_BOUNDARY_HH__ */

