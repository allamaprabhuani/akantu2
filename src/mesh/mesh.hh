/**
 * @file   mesh.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  the class representing the meshes
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
#ifndef __AKANTU_MESH_HH__
#define __AKANTU_MESH_HH__

/* -------------------------------------------------------------------------- */

#include "aka_config.hh"
#include "aka_common.hh"
#include "aka_memory.hh"
#include "aka_array.hh"
#include "element_class.hh"
#include "element_type_map.hh"
#include "aka_event_handler_manager.hh"
#include "group_manager.hh"
#include "element.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */
#include "mesh_data.hh"
#include "dumpable.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Mesh modifications events                                                  */
/* -------------------------------------------------------------------------- */
#include "mesh_events.hh"


/* -------------------------------------------------------------------------- */
/* Mesh                                                                       */
/* -------------------------------------------------------------------------- */

/**
 * @class  Mesh this  contain the  coordinates of  the nodes  in  the Mesh.nodes
 * Array,  and the  connectivity. The  connectivity  are stored  in by  element
 * types.
 *
 * To  know  all  the  element  types   present  in  a  mesh  you  can  get  the
 * Mesh::ConnectivityTypeList
 *
 * In order to loop on all element you have to loop on all types like this :
 * @code{.cpp}
 Mesh::type_iterator it = mesh.firstType(dim, ghost_type);
 Mesh::type_iterator end = mesh.lastType(dim, ghost_type);
 for(; it != end; ++it) {
   UInt nb_element  = mesh.getNbElement(*it);
   const Array<UInt> & conn = mesh.getConnectivity(*it);
   for(UInt e = 0; e < nb_element; ++e) {
     ...
   }
 }
 @endcode
*/
class Mesh : protected Memory,
             public EventHandlerManager<MeshEventHandler>,
             public GroupManager,
	     public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// constructor that create nodes coordinates array
  Mesh(UInt spatial_dimension,
       const ID & id = "mesh",
       const MemoryID & memory_id = 0);

  /// constructor that use an existing nodes coordinates array, by knowing its ID
  Mesh(UInt spatial_dimension,
       const ID & nodes_id,
       const ID & id,
       const MemoryID & memory_id = 0);

  /**
   * constructor that use an existing nodes coordinates
   * array, by getting the vector of coordinates
   */
  Mesh(UInt spatial_dimension,
       Array<Real> & nodes,
       const ID & id = "mesh",
       const MemoryID & memory_id = 0);


  virtual ~Mesh();

  /// @typedef ConnectivityTypeList list of the types present in a Mesh
  typedef std::set<ElementType> ConnectivityTypeList;

  /// read the mesh from a file
  void read (const std::string & filename, const MeshIOType & mesh_io_type = _miot_auto);
  /// write the mesh to a file
  void write(const std::string & filename, const MeshIOType & mesh_io_type = _miot_auto);

private:
  /// initialize the connectivity to NULL and other stuff
  void init();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// function that computes the bounding box (fills xmin, xmax)
  void computeBoundingBox();

#ifdef AKANTU_CORE_CXX11
  /// translate the mesh by the given amount in x[, y[, z]] directions
  template <typename... Args> void translate(Args... params);
#endif

  /// init a by-element-type real vector with provided ids
  template<typename T>
  void initElementTypeMapArray(ElementTypeMapArray<T> & v,
                              UInt nb_component,
                              UInt spatial_dimension,
                              const bool & flag_nb_node_per_elem_multiply = false,
                              ElementKind element_kind = _ek_regular,
                              bool size_to_nb_element = false) const; /// @todo: think about nicer way to do it
  template<typename T>
  void initElementTypeMapArray(ElementTypeMapArray<T> & v,
                              UInt nb_component,
                              UInt spatial_dimension,
                              GhostType ghost_type,
                              const bool & flag_nb_node_per_elem_multiply = false,
                              ElementKind element_kind = _ek_regular,
                              bool size_to_nb_element = false) const; /// @todo: think about nicer way to do it

  template<typename T>
  void initElementTypeMapArray(ElementTypeMapArray<T> & v,
			       UInt nb_component,
			       UInt spatial_dimension,
			       GhostType ghost_type,
			       const T & default_value,
			       const bool & flag_nb_node_per_elem_multiply = false,
			       ElementKind element_kind = _ek_regular,
			       bool size_to_nb_element = false) const; /// @todo: think about nicer way to do it

  /// extract coordinates of nodes from an element
  template<typename T>
  inline void extractNodalValuesFromElement(const Array<T> & nodal_values,
                                            T * elemental_values,
                                            UInt * connectivity,
                                            UInt n_nodes,
                                            UInt nb_degree_of_freedom) const;

  // /// extract coordinates of nodes from a reversed element
  // inline void extractNodalCoordinatesFromPBCElement(Real * local_coords,
  //                                                   UInt * connectivity,
  //                                                   UInt n_nodes);

  /// convert a element to a linearized element
  inline UInt elementToLinearized(const Element & elem) const;

  /// convert a linearized element to an element
  inline Element linearizedToElement (UInt linearized_element) const;

  /// update the types offsets array for the conversions
  inline void updateTypesOffsets(const GhostType & ghost_type);

  /// add a Array of connectivity for the type <type>.
  inline void addConnectivityType(const ElementType & type,
                                  const GhostType & ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  template <class Event>
  inline void sendEvent(Event & event) {
    //    if(event.getList().getSize() != 0)
    EventHandlerManager<MeshEventHandler>::sendEvent<Event>(event);
  }

  /* ------------------------------------------------------------------------ */
  template<typename T>
  inline void removeNodesFromArray(Array<T> & vect, const Array<UInt> & new_numbering);

  /// initialize normals
  void initNormals();

  /// init facets' mesh
  Mesh & initMeshFacets(const ID & id = "mesh_facets");

  /// define parent mesh
  void defineMeshParent(const Mesh & mesh);

  /// get global connectivity array
  void getGlobalConnectivity(ElementTypeMapArray<UInt> & global_connectivity,
			     UInt dimension,
			     GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the id of the mesh
  AKANTU_GET_MACRO(ID, Memory::id, const ID &);

  /// get the spatial dimension of the mesh = number of component of the coordinates
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the nodes Array aka coordinates
  AKANTU_GET_MACRO(Nodes, *nodes, const Array<Real> &);
  AKANTU_GET_MACRO_NOT_CONST(Nodes, *nodes, Array<Real> &);

  /// get the normals for the elements
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Normals, normals, Real);

  /// get the number of nodes
  AKANTU_GET_MACRO(NbNodes, nodes->getSize(), UInt);

  /// get the Array of global ids of the nodes (only used in parallel)
  AKANTU_GET_MACRO(GlobalNodesIds, *nodes_global_ids, const Array<UInt> &);
  AKANTU_GET_MACRO_NOT_CONST(GlobalNodesIds, *nodes_global_ids, Array<UInt> &);

  /// get the global id of a node
  inline UInt getNodeGlobalId(UInt local_id) const;

  /// get the global number of nodes
  inline UInt getNbGlobalNodes() const;

  /// get the nodes type Array
  AKANTU_GET_MACRO(NodesType, nodes_type, const Array<Int> &);

protected:
  AKANTU_GET_MACRO_NOT_CONST(NodesType, nodes_type, Array<Int> &);
public:
  inline Int getNodeType(UInt local_id) const;

  /// say if a node is a pure ghost node
  inline bool isPureGhostNode(UInt n) const;

  /// say if a node is pur local or master node
  inline bool isLocalOrMasterNode(UInt n) const;

  inline bool isLocalNode(UInt n) const;
  inline bool isMasterNode(UInt n) const;
  inline bool isSlaveNode(UInt n) const;

  AKANTU_GET_MACRO(LowerBounds, lower_bounds, const Vector<Real> &);
  AKANTU_GET_MACRO(UpperBounds, upper_bounds, const Vector<Real> &);
  AKANTU_GET_MACRO(LocalLowerBounds, local_lower_bounds, const Vector<Real> &);
  AKANTU_GET_MACRO(LocalUpperBounds, local_upper_bounds, const Vector<Real> &);

  /// get the connectivity Array for a given type
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Connectivity, connectivities, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Connectivity, connectivities, UInt);
  AKANTU_GET_MACRO(Connectivities, connectivities, const ElementTypeMapArray<UInt> &);

  /// get the number of element of a type in the mesh
  inline UInt getNbElement(const ElementType & type, const GhostType & ghost_type = _not_ghost) const;

  /// get the number of element for a given ghost_type and a given dimension
  inline UInt getNbElement(const UInt spatial_dimension = _all_dimensions,
			   const GhostType & ghost_type = _not_ghost,
			   const ElementKind & kind = _ek_not_defined) const;

  /// get the connectivity list either for the elements or the ghost elements
  inline const ConnectivityTypeList & getConnectivityTypeList(const GhostType & ghost_type = _not_ghost) const;

  /// compute the barycenter of a given element
  inline void getBarycenter(UInt element, const ElementType & type, Real * barycenter,
                            GhostType ghost_type = _not_ghost) const;
  inline void getBarycenter(const Element & element, Vector<Real> & barycenter) const;

  /// get the element connected to a subelement
  const Array< std::vector<Element> > & getElementToSubelement(const ElementType & el_type,
                                                               const GhostType & ghost_type = _not_ghost) const;
  /// get the element connected to a subelement
  Array< std::vector<Element> > & getElementToSubelement(const ElementType & el_type,
                                                         const GhostType & ghost_type = _not_ghost);

  /// get the subelement connected to an element
  const Array<Element> & getSubelementToElement(const ElementType & el_type,
                                                const GhostType & ghost_type = _not_ghost) const;
  /// get the subelement connected to an element
  Array<Element> & getSubelementToElement(const ElementType & el_type,
                                          const GhostType & ghost_type = _not_ghost);

  /// get a name field associated to the mesh
  template<typename T>
  inline const Array<T> & getData(const std::string & data_name,
				  const ElementType & el_type,
                                  const GhostType & ghost_type = _not_ghost) const;

  /// get a name field associated to the mesh
  template<typename T>
  inline Array<T> & getData(const std::string & data_name,
			    const ElementType & el_type,
                            const GhostType & ghost_type = _not_ghost);

  /// register a new ElementalTypeMap in the MeshData
  template<typename T>
  inline ElementTypeMapArray<T> & registerData(const std::string & data_name);

  /// get a name field associated to the mesh
  template<typename T>
  inline const ElementTypeMapArray<T> & getData(const std::string & data_name) const;

  /// get a name field associated to the mesh
  template<typename T>
  inline ElementTypeMapArray<T> & getData(const std::string & data_name);


  template <typename T>
  ElementTypeMap<UInt> getNbDataPerElem(ElementTypeMapArray<T> & array,
					const ElementKind & element_kind);

  template <typename T>
  dumper::Field * createFieldFromAttachedData(const std::string & field_id,
					      const std::string & group_name,
					      const ElementKind & element_kind);

  /// templated getter returning the pointer to data in MeshData (modifiable)
  template<typename T>
  inline Array<T> * getDataPointer(const std::string & data_name,
				   const ElementType & el_type,
                                   const GhostType & ghost_type = _not_ghost,
                                   UInt nb_component = 1,
				   bool size_to_nb_element = true,
				   bool resize_with_parent = false);

  /// Facets mesh accessor
  AKANTU_GET_MACRO(MeshFacets, *mesh_facets, const Mesh &);
  AKANTU_GET_MACRO_NOT_CONST(MeshFacets, *mesh_facets, Mesh &);

  /// Parent mesh accessor
  AKANTU_GET_MACRO(MeshParent, *mesh_parent, const Mesh &);

  inline bool isMeshFacets() const { return this->is_mesh_facets; }

  /// defines is the mesh is distributed or not
  inline bool isDistributed() const { return this->is_distributed; }

#ifndef SWIG
  /// return the dumper from a group and and a dumper name
  DumperIOHelper & getGroupDumper(const std::string & dumper_name, 
				  const std::string & group_name);
#endif
  /* ------------------------------------------------------------------------ */
  /* Wrappers on ElementClass functions                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of nodes per element for a given element type
  static inline UInt getNbNodesPerElement(const ElementType & type);

  /// get the number of nodes per element for a given element type considered as
  /// a first order element
  static inline ElementType getP1ElementType(const ElementType & type);

  /// get the kind of the element type
  static inline ElementKind getKind(const ElementType & type);

  /// get spatial dimension of a type of element
  static inline UInt getSpatialDimension(const ElementType & type);

  /// get number of facets of a given element type
  static inline UInt getNbFacetsPerElement(const ElementType & type);

  /// get number of facets of a given element type
  static inline UInt getNbFacetsPerElement(const ElementType & type, UInt t);

  /// get local connectivity of a facet for a given facet type
  static inline MatrixProxy<UInt> getFacetLocalConnectivity(const ElementType & type, UInt t = 0);

  /// get connectivity of facets for a given element
  inline Matrix<UInt> getFacetConnectivity(const Element & element, UInt t = 0) const;

  /// get the number of type of the surface element associated to a given element type
  static inline UInt getNbFacetTypes(const ElementType & type, UInt t = 0);

  /// get the type of the surface element associated to a given element
  static inline ElementType getFacetType(const ElementType & type, UInt t = 0);

  /// get all the type of the surface element associated to a given element
  static inline VectorProxy<ElementType> getAllFacetTypes(const ElementType & type);

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  typedef ElementTypeMapArray<UInt, ElementType>::type_iterator type_iterator;

  inline type_iterator firstType(UInt dim = _all_dimensions,
                                 GhostType ghost_type = _not_ghost,
                                 ElementKind kind = _ek_regular) const {
    return connectivities.firstType(dim, ghost_type, kind);
  }

  inline type_iterator lastType(UInt dim = _all_dimensions,
                                GhostType ghost_type = _not_ghost,
                                ElementKind kind = _ek_regular) const {
    return connectivities.lastType(dim, ghost_type, kind);
  }

  /* ------------------------------------------------------------------------ */
  /* Private methods for friends                                              */
  /* ------------------------------------------------------------------------ */
private:
  friend class MeshAccessor;

  friend class MeshIOMSH;
  friend class MeshIOMSHStruct;
  friend class MeshIODiana;
  friend class MeshUtils;
  friend class DistributedSynchronizer;
  template<class T> friend class SpatialGrid;

#if defined(AKANTU_COHESIVE_ELEMENT)
  friend class CohesiveElementInserter;
#endif

#if defined(AKANTU_IGFEM)
  template<UInt dim> friend class MeshIgfemSphericalGrowingGel;
#endif

  AKANTU_GET_MACRO(NodesPointer, nodes, Array<Real> *);

  /// get a pointer to the nodes_global_ids Array<UInt> and create it if necessary
  inline Array<UInt> * getNodesGlobalIdsPointer();

  /// get a pointer to the nodes_type Array<Int> and create it if necessary
  inline Array<Int> * getNodesTypePointer();

  /// get a pointer to the connectivity Array for the given type and create it if necessary
  inline Array<UInt> * getConnectivityPointer(const ElementType & type,
                                              const GhostType & ghost_type = _not_ghost);



  /// get a pointer to the element_to_subelement Array for the given type and create it if necessary
  inline Array< std::vector<Element> > * getElementToSubelementPointer(const ElementType & type,
                                                                       const GhostType & ghost_type = _not_ghost);

  /// get a pointer to the subelement_to_element Array for the given type and create it if necessary
  inline Array<Element > * getSubelementToElementPointer(const ElementType & type,
                                                         const GhostType & ghost_type = _not_ghost);

  AKANTU_GET_MACRO_NOT_CONST(MeshData, mesh_data, MeshData &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// array of the nodes coordinates
  Array<Real> * nodes;

  /// global node ids
  Array<UInt> * nodes_global_ids;

  /// node type,  -3 pure ghost, -2  master for the  node, -1 normal node,  i in
  /// [0-N] slave node and master is proc i
  Array<Int> nodes_type;

  /// global number of nodes;
  UInt nb_global_nodes;

  /// boolean to know if the nodes have to be deleted with the mesh or not
  bool created_nodes;

  /// all class of elements present in this mesh (for heterogenous meshes)
  ElementTypeMapArray<UInt> connectivities;

  /// map to normals for all class of elements present in this mesh
  ElementTypeMapArray<Real> normals;

  /// list of all existing types in the mesh
  ConnectivityTypeList type_set;

  /// the spatial dimension of this mesh
  UInt spatial_dimension;

  /// types offsets
  Array<UInt> types_offsets;

  /// list of all existing types in the mesh
  ConnectivityTypeList ghost_type_set;
  /// ghost types offsets
  Array<UInt> ghost_types_offsets;

  /// min of coordinates
  Vector<Real> lower_bounds;
  /// max of coordinates
  Vector<Real> upper_bounds;
  /// size covered by the mesh on each direction
  Vector<Real> size;

  /// local min of coordinates
  Vector<Real> local_lower_bounds;
  /// local max of coordinates
  Vector<Real> local_upper_bounds;

  /// Extra data loaded from the mesh file
  MeshData mesh_data;

  /// facets' mesh
  Mesh * mesh_facets;

  /// parent mesh (this is set for mesh_facets meshes)
  const Mesh * mesh_parent;

  /// defines if current mesh is mesh_facets or not
  bool is_mesh_facets;

  /// defines if the mesh is centralized or distributed
  bool is_distributed;
};


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Element & _this)
{
  _this.printself(stream);
  return stream;
}


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Mesh & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__


/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "mesh_inline_impl.cc"
#include "element_type_map_tmpl.hh"
//#include "group_manager_inline_impl.cc"
//#include "element_group_inline_impl.cc"

#endif /* __AKANTU_MESH_HH__ */
