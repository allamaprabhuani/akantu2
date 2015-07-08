/**
 * @file   group_manager.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  Stores information relevent to the notion of element and nodes groups.
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

#ifndef __AKANTU_GROUP_MANAGER_HH__
#define __AKANTU_GROUP_MANAGER_HH__

#include <set>
#include "aka_common.hh"
#include "element_type_map.hh"
//#include "dumpable.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class FEM;
class ElementGroup;
class NodeGroup;
class Mesh;
class Element;
class DistributedSynchronizer;
template<bool> class CommunicationBufferTemplated;

namespace dumper {
  class Field;
}
/* -------------------------------------------------------------------------- */



class GroupManager {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::map<std::string, ElementGroup *> ElementGroups;
  typedef std::map<std::string, NodeGroup *> NodeGroups;

public:
  typedef std::set<ElementType> GroupManagerTypeSet;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  GroupManager(const Mesh & mesh,
	       const ID & id = "group_manager",
	       const MemoryID & memory_id = 0);
  virtual ~GroupManager();

  /* ------------------------------------------------------------------------ */
  /* Groups iterators                                                         */
  /* ------------------------------------------------------------------------ */
public:
  typedef NodeGroups::iterator node_group_iterator;
  typedef ElementGroups::iterator element_group_iterator;

  typedef NodeGroups::const_iterator const_node_group_iterator;
  typedef ElementGroups::const_iterator const_element_group_iterator;

#ifndef SWIG
#define AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION(group_type,       \
						      function,         \
						      param_in,         \
						      param_out)        \
  inline BOOST_PP_CAT(BOOST_PP_CAT(const_, group_type), _iterator)      \
    BOOST_PP_CAT(BOOST_PP_CAT(group_type, _), function)(param_in) const { \
    return BOOST_PP_CAT(group_type, s).function(param_out);             \
  };                                                                    \
									\
  inline BOOST_PP_CAT(group_type, _iterator)                            \
    BOOST_PP_CAT(BOOST_PP_CAT(group_type, _), function)(param_in) {     \
    return BOOST_PP_CAT(group_type, s).function(param_out);             \
  }

#define AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(group_type, function) \
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION(group_type, function, BOOST_PP_EMPTY(), BOOST_PP_EMPTY())

  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(node_group, begin);
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(node_group, end  );
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(element_group, begin);
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(element_group, end  );
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION(element_group, find, const std::string & name, name);
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION(node_group, find, const std::string & name, name);
#endif

  /* ------------------------------------------------------------------------ */
  /* Clustering filter                                                        */
  /* ------------------------------------------------------------------------ */
public:

  class ClusteringFilter {
  public:
    virtual bool operator() (const Element &) const {
      return true;
    }
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// create an empty node group
  NodeGroup & createNodeGroup(const std::string & group_name,
			      bool replace_group = false);
  
  /// create a node group from another node group but filtered
  template <typename T>
  NodeGroup & createFilteredNodeGroup(const std::string & group_name,
				      const NodeGroup & node_group,
				      T & filter);

  /// destroy a node group
  void destroyNodeGroup(const std::string & group_name);


  /// create an element group and the associated node group
  ElementGroup & createElementGroup(const std::string & group_name,
				    UInt dimension = _all_dimensions,
				    bool replace_group = false);

  /// create an element group from another element group but filtered
  template <typename T>
  ElementGroup & createFilteredElementGroup(const std::string & group_name,
					    UInt dimension,
					    const NodeGroup & node_group,
					    T & filter);

  /// destroy an element group and the associated node group
  void destroyElementGroup(const std::string & group_name,
			   bool destroy_node_group = false);

  /// destroy all element groups and the associated node groups
  void destroyAllElementGroups(bool destroy_node_groups = false);

  /// create a element group using an existing node group
  ElementGroup & createElementGroup(const std::string & group_name,
				    UInt dimension,
				    NodeGroup & node_group);

  /// create groups based on values stored in a given mesh data
  template <typename T>
  void createGroupsFromMeshData(const std::string & dataset_name);

  /// create boundaries group by a clustering algorithm \todo extend to parallel
  UInt createBoundaryGroupFromGeometry();

  /// create element clusters for a given dimension
  UInt createClusters(UInt element_dimension,
		      std::string cluster_name_prefix = "cluster",
		      const ClusteringFilter & filter = ClusteringFilter(),
		      DistributedSynchronizer * distributed_synchronizer = NULL,
		      Mesh * mesh_facets = NULL);

  /// Create an ElementGroup based on a NodeGroup
  void createElementGroupFromNodeGroup(const std::string & name,
				       const std::string & node_group,
				       UInt dimension = _all_dimensions);

  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// this function insure that the group names are present on all processors
  /// /!\ it is a SMP call
  void synchronizeGroupNames();

  /// register an elemental field to the given group name (overloading for ElementalPartionField)
#ifndef SWIG
  template <typename T, template <bool> class dump_type>
  inline dumper::Field * createElementalField(const ElementTypeMapArray<T> & field,
					      const std::string & group_name,
					      UInt spatial_dimension,
					      const ElementKind & kind,
					      ElementTypeMap<UInt> nb_data_per_elem = ElementTypeMap<UInt>());

  /// register an elemental field to the given group name (overloading for ElementalField)
  template <typename T, template <class> class ret_type, template <class,template <class> class,bool> class dump_type>
  inline dumper::Field * createElementalField(const ElementTypeMapArray<T> & field,
					      const std::string & group_name,
					      UInt spatial_dimension,
					      const ElementKind & kind,
					      ElementTypeMap<UInt> nb_data_per_elem = ElementTypeMap<UInt>());

  /// register an elemental field to the given group name (overloading for MaterialInternalField)
  template <typename T,
	    /// type of InternalMaterialField
	    template<typename T, bool filtered> class dump_type>
  inline dumper::Field * createElementalField(const ElementTypeMapArray<T> & field,
					      const std::string & group_name,
					      UInt spatial_dimension,
					      const ElementKind & kind,
					      ElementTypeMap<UInt> nb_data_per_elem);

  template <typename type, bool flag, template<class,bool> class ftype>
  inline dumper::Field * createNodalField(const ftype<type,flag> * field,
					  const std::string & group_name,
					  UInt padding_size = 0);

  template <typename type, bool flag, template<class,bool> class ftype>
  inline dumper::Field * createStridedNodalField(const ftype<type,flag> * field,
						 const std::string & group_name,
						 UInt size, UInt stride,
						 UInt padding_size);

protected:
  /// fill a buffer with all the group names
  void fillBufferWithGroupNames(CommunicationBufferTemplated<false> & comm_buffer) const;

  /// take a buffer and create the missing groups localy
  void checkAndAddGroups(CommunicationBufferTemplated<true> & buffer);



  /// register an elemental field to the given group name
  template <class dump_type,typename field_type>
  inline dumper::Field * createElementalField(const field_type & field,
				       const std::string & group_name,
				       UInt spatial_dimension,
				       const ElementKind & kind,
				       ElementTypeMap<UInt> nb_data_per_elem);


  /// register an elemental field to the given group name
  template <class dump_type,typename field_type>
  inline dumper::Field * createElementalFilteredField(const field_type & field,
					       const std::string & group_name,
					       UInt spatial_dimension,
					       const ElementKind & kind,
					       ElementTypeMap<UInt> nb_data_per_elem);
#endif

  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  const ElementGroup & getElementGroup(const std::string & name) const;
  const NodeGroup    & getNodeGroup(const std::string & name) const;

  ElementGroup & getElementGroup(const std::string & name);
  NodeGroup    & getNodeGroup(const std::string & name);


  UInt getNbElementGroups(UInt dimension = _all_dimensions) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id to create element and node groups
  ID id;
  /// memory_id to create element and node groups
  MemoryID memory_id;

  /// list of the node groups managed
  NodeGroups node_groups;

  /// list of the element groups managed
  ElementGroups element_groups;

  /// Mesh to which the element belongs
  const Mesh & mesh;
};

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const GroupManager & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_GROUP_MANAGER_HH__ */
