/**
 * @file   mesh.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 *
 * @date   Fri Jun 18 11:47:19 2010
 *
 * @brief  the class representing the meshes
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
#ifndef __AKANTU_MESH_HH__
#define __AKANTU_MESH_HH__

/* -------------------------------------------------------------------------- */

#include "aka_config.hh"
#include "aka_common.hh"
#include "aka_memory.hh"
#include "aka_vector.hh"
#include "element_class.hh"
#include "by_element_type.hh"
#include "aka_event_handler.hh"
#include "boundary.hh"

/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class SubBoundary;

/* -------------------------------------------------------------------------- */
/* Element                                                                    */
/* -------------------------------------------------------------------------- */
class Element;
extern const Element ElementNull;

class Element {
public:
  Element(ElementType type = _not_defined, UInt element = 0,
          GhostType ghost_type = _not_ghost, ElementKind kind = _ek_regular) :
    type(type), element(element),
    ghost_type(ghost_type), kind(kind) {};

  Element(const Element & element) {
    this->type    = element.type;
    this->element = element.element;
    this->ghost_type = element.ghost_type;
    this->kind = element.kind;
  }

  inline bool operator==(const Element & elem) const {
    return ((element == elem.element)
            && (type == elem.type)
            && (ghost_type == elem.ghost_type)
            && (kind == elem.kind));
  }

  inline bool operator!=(const Element & elem) const {
    return ((element != elem.element)
            || (type != elem.type)
            || (ghost_type != elem.ghost_type)
            || (kind != elem.kind));
  }

  bool operator<(const Element& rhs) const {
    bool res = (rhs == ElementNull) || ((this->kind < rhs.kind) ||
                                        ((this->kind == rhs.kind) &&
                                         ((this->ghost_type < rhs.ghost_type) ||
                                          ((this->ghost_type == rhs.ghost_type) &&
                                           ((this->type < rhs.type) ||
                                            ((this->type == rhs.type) &&
                                             (this->element < rhs.element)))))));
    return res;
  }

  virtual ~Element() {};

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
public:
  ElementType type;
  UInt element;
  GhostType ghost_type;
  ElementKind kind;
};

struct CompElementLess {
  bool operator() (const Element& lhs, const Element& rhs) const {
    return lhs < rhs;
  }
};

__END_AKANTU__
#include "mesh_data.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Mesh modifications events                                                  */
/* -------------------------------------------------------------------------- */
template<class Entity>
class MeshEvent {
public:
  virtual ~MeshEvent() {}
  const Array<Entity> & getList() const { return list; }
  Array<Entity> & getList() { return list; }
protected:
  Array<Entity> list;
};

class Mesh;

class NewNodesEvent : public MeshEvent<UInt> {
public:
  virtual ~NewNodesEvent() {};
};
class RemovedNodesEvent : public MeshEvent<UInt> {
public:
  virtual ~RemovedNodesEvent() {};
  inline RemovedNodesEvent(const Mesh & mesh);
  AKANTU_GET_MACRO_NOT_CONST(NewNumbering, new_numbering, Array<UInt> &);
  AKANTU_GET_MACRO(NewNumbering, new_numbering, const Array<UInt> &);
private:
  Array<UInt> new_numbering;
};

class NewElementsEvent : public MeshEvent<Element> {
public:
  virtual ~NewElementsEvent() {};
};
class RemovedElementsEvent : public MeshEvent<Element> {
public:
  virtual ~RemovedElementsEvent() {};
  inline RemovedElementsEvent(const Mesh & mesh);
  AKANTU_GET_MACRO(NewNumbering, new_numbering, const ByElementTypeUInt &);
  AKANTU_GET_MACRO_NOT_CONST(NewNumbering, new_numbering, ByElementTypeUInt &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(NewNumbering, new_numbering, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(NewNumbering, new_numbering, UInt);
protected:
  ByElementTypeUInt new_numbering;
};

/* -------------------------------------------------------------------------- */

class MeshEventHandler {
public:
  virtual ~MeshEventHandler() {};
  /* ------------------------------------------------------------------------ */
  /* Internal code                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  inline void sendEvent(const NewNodesEvent & event)     { onNodesAdded  (event.getList(),
                                                                          event); }
  inline void sendEvent(const RemovedNodesEvent & event) { onNodesRemoved(event.getList(),
                                                                          event.getNewNumbering(),
                                                                          event); }

  inline void sendEvent(const NewElementsEvent & event)     { onElementsAdded  (event.getList(),
                                                                                event); }
  inline void sendEvent(const RemovedElementsEvent & event) { onElementsRemoved(event.getList(),
                                                                                event.getNewNumbering(),
                                                                                event); }

  template<class EventHandler>
  friend class EventHandlerManager;

  /* ------------------------------------------------------------------------ */
  /* Interface                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onNodesAdded  (__attribute__((unused)) const Array<UInt> & nodes_list,
                              __attribute__((unused)) const NewNodesEvent & event) {  }
  virtual void onNodesRemoved(__attribute__((unused)) const Array<UInt> & nodes_list,
                              __attribute__((unused)) const Array<UInt> & new_numbering,
                              __attribute__((unused)) const RemovedNodesEvent & event) {  }

  virtual void onElementsAdded  (__attribute__((unused)) const Array<Element> & elements_list,
                                 __attribute__((unused)) const NewElementsEvent & event) { }
  virtual void onElementsRemoved(__attribute__((unused)) const Array<Element> & elements_list,
                                 __attribute__((unused)) const ByElementTypeUInt & new_numbering,
                                 __attribute__((unused)) const RemovedElementsEvent & event) { }
};

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
 * @code
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
class Mesh : protected Memory, public EventHandlerManager<MeshEventHandler> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// constructor that create nodes coordinates array
  Mesh(UInt spatial_dimension,
       const ID id = "mesh",
       const MemoryID & memory_id = 0);

  /// constructor that use an existing nodes coordinates array, by knowing its ID
  Mesh(UInt spatial_dimension,
       const ID & nodes_id,
       const ID & id,
       const MemoryID & memory_id);

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
  template <typename... Args>
  void translate(Args... params) {
    // check that the number of parameters corresponds to the dimension
    AKANTU_DEBUG_ASSERT(sizeof...(Args) <= spatial_dimension , "Number of arguments greater than dimension.");

    // unpack parameters
    Real s[] = { params... };

    Array<Real>& nodes = getNodes();
    for (UInt i = 0; i < nodes.getSize(); ++i)
      for (UInt k = 0; k < sizeof...(Args); ++k)
        nodes(i, k) += s[k];
  }
#endif

  /// init a by-element-type real vector with provided ids
  template<typename T>
  void initByElementTypeArray(ByElementTypeArray<T> & v,
                              UInt nb_component,
                              UInt spatial_dimension,
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

  /// extract coordinates of nodes from a reversed element
  inline void extractNodalCoordinatesFromPBCElement(Real * local_coords,
                                                    UInt * connectivity,
                                                    UInt n_nodes);

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

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(ID, id, const ID &);

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

  /// get the global id of a node
  inline UInt getNodeGlobalId(UInt local_id) const;

  /// get the global number of nodes
  inline UInt getNbGlobalNodes() const;

  /// get the nodes type Array
  AKANTU_GET_MACRO(NodesType, *nodes_type, const Array<Int> &);
  inline Int getNodeType(UInt local_id) const;

  /// say if a node is a pure ghost node
  inline bool isPureGhostNode(UInt n) const;

  /// say if a node is pur local or master node
  inline bool isLocalOrMasterNode(UInt n) const;

  inline bool isLocalNode(UInt n) const;
  inline bool isMasterNode(UInt n) const;
  inline bool isSlaveNode(UInt n) const;

  AKANTU_GET_MACRO(XMin, lower_bounds[0], Real);
  AKANTU_GET_MACRO(YMin, lower_bounds[1], Real);
  AKANTU_GET_MACRO(ZMin, lower_bounds[2], Real);
  AKANTU_GET_MACRO(XMax, upper_bounds[0], Real);
  AKANTU_GET_MACRO(YMax, upper_bounds[1], Real);
  AKANTU_GET_MACRO(ZMax, upper_bounds[2], Real);

  inline void getLowerBounds(Real * lower) const;
  inline void getUpperBounds(Real * upper) const;

  inline void getLocalLowerBounds(Real * lower) const;
  inline void getLocalUpperBounds(Real * upper) const;

  /// get the connectivity Array for a given type
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Connectivity, connectivities, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Connectivity, connectivities, UInt);
  AKANTU_GET_MACRO(Connectivities, connectivities, const ByElementTypeArray<UInt> &);

  /// get the number of element of a type in the mesh
  inline UInt getNbElement(const ElementType & type, const GhostType & ghost_type = _not_ghost) const;

  /// get the number of element for a given ghost_type and a given dimension
  inline UInt getNbElement(const UInt spatial_dimension = _all_dimensions, const GhostType & ghost_type = _not_ghost) const;

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


  inline UInt getNbBoundaries() const;
  inline const SubBoundary & getSubBoundary(const std::string & name) const;
  inline SubBoundary & getSubBoundary(const std::string & name);

  // // XXX TODO FIXME old prototype
  // inline const Array<UInt> & getUIntData(const ElementType & el_type,
  //                                        const std::string & data_name,
  //                                        const GhostType & ghost_type = _not_ghost) const;

  AKANTU_GET_MACRO(Boundary, boundaries, const Boundary &);
  AKANTU_GET_MACRO_NOT_CONST(Boundary, boundaries, Boundary &);

  /// get a name field associated to the mesh
  template<typename T>
  inline const Array<T> & getData(const ElementType & el_type,
                                  const std::string & data_name,
                                  const GhostType & ghost_type = _not_ghost) const;

  /// get a name field associated to the mesh
  template<typename T>
  inline Array<T> & getData(const ElementType & el_type,
                            const std::string & data_name,
                            const GhostType & ghost_type = _not_ghost);

  /// get a name field associated to the mesh
  template<typename T>
  inline const ByElementTypeArray<T> & getData(const std::string & data_name) const;

  /// get a name field associated to the mesh
  template<typename T>
  inline ByElementTypeArray<T> & getData(const std::string & data_name);

  /// templated getter returning the pointer to data in MeshData (modifiable)
  template<typename T>
  inline Array<T> * getDataPointer(const ElementType & el_type,
                                   const std::string & data_name,
                                   const GhostType & ghost_type = _not_ghost,
                                   UInt nb_component = 1);

  /// Facets mesh accessor
  AKANTU_GET_MACRO(MeshFacets, *mesh_facets, const Mesh &);
  AKANTU_GET_MACRO_NOT_CONST(MeshFacets, *mesh_facets, Mesh &);

  /// Parent mesh accessor
  AKANTU_GET_MACRO(MeshParent, *mesh_parent, const Mesh &);
  AKANTU_GET_MACRO_NOT_CONST(MeshParent, *mesh_parent, Mesh &);

  inline bool isMeshFacets() const {return is_mesh_facets;}

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

  /// get local connectivity of a facet for a given facet type
  static inline Matrix<UInt> getFacetLocalConnectivity(const ElementType & type);

  /// get connectivity of facets for a given element
  inline Matrix<UInt> getFacetConnectivity(UInt element, const ElementType & type, const GhostType & ghost_type) const;

  /// get the type of the surface element associated to a given element
  static inline ElementType getFacetType(const ElementType & type);

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  typedef ByElementTypeArray<UInt, ElementType>::type_iterator type_iterator;

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
  friend class MeshIOMSH;
  friend class MeshIOMSHStruct;
  friend class MeshIODiana;
  friend class MeshUtils;
  friend class DistributedSynchronizer;
  template<class T> friend class SpatialGrid;


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
  /// id of the mesh
  ID id;

  /// array of the nodes coordinates
  Array<Real> * nodes;

  /// global node ids
  Array<UInt> * nodes_global_ids;

  /// node type,  -3 pure ghost, -2  master for the  node, -1 normal node,  i in
  /// [0-N] slave node and master is proc i
  Array<Int> * nodes_type;

  /// global number of nodes;
  UInt nb_global_nodes;

  /// boolean to know if the nodes have to be deleted with the mesh or not
  bool created_nodes;

  /// all class of elements present in this mesh (for heterogenous meshes)
  ByElementTypeUInt connectivities;

  /// map to normals for all class of elements present in this mesh
  ByElementTypeReal normals;

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
  Real lower_bounds[3];
  /// max of coordinates
  Real upper_bounds[3];
  /// size covered by the mesh on each direction
  Real size[3];

  /// local min of coordinates
  Real local_lower_bounds[3];
  /// local max of coordinates
  Real local_upper_bounds[3];

  /// Extra data loaded from the mesh file
  MeshData mesh_data;

  /// List of boundaries either read from the mesh file or created directly
  Boundary boundaries;

  /// facets' mesh
  Mesh * mesh_facets;

  /// parent mesh (this is set for mesh_facets meshes)
  Mesh * mesh_parent;

  /// defines if current mesh is mesh_facets or not
  bool is_mesh_facets;

};

/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Element & _this)
{
  _this.printself(stream);
  return stream;
}

#include "mesh_inline_impl.cc"

#include "by_element_type_tmpl.hh"


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Mesh & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__


#endif /* __AKANTU_MESH_HH__ */
