/**
 * @file   mesh.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 11:53:53 2010
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
#include "aka_common.hh"
#include "aka_memory.hh"
#include "aka_vector.hh"
#include "element_class.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Element                                                                    */
/* -------------------------------------------------------------------------- */
class Element {
public:
  Element(ElementType type = _not_defined, UInt element = 0, GhostType ghost_type = _not_ghost) :
    type(type), element(element), ghost_type(ghost_type) {};

  Element(const Element & element) {
    this->type    = element.type;
    this->element = element.element;
    this->ghost_type = element.ghost_type;
  }

  virtual ~Element() {};

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
public:
  ElementType type;
  UInt element;
  GhostType ghost_type;
};

struct CompElementLess {
  bool operator() (const Element& lhs, const Element& rhs) const {
    bool res = ((lhs.ghost_type < rhs.ghost_type) ||
                ((lhs.ghost_type == rhs.ghost_type) &&
                 ((lhs.type < rhs.type) ||
                  ((lhs.type == rhs.type) &&
                   (lhs.element < rhs.element)))));
    return res;
  }
};

extern const Element ElementNull;

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

  inline static std::string printType(const ElementType & type, const GhostType & ghost_type);

  inline bool exists(ElementType type, GhostType ghost_type = _not_ghost) const;

  inline const Stored & operator()(const ElementType & type,
				   const GhostType & ghost_type = _not_ghost) const;
  inline Stored & operator()(const ElementType & type,
			     const GhostType & ghost_type = _not_ghost);

  inline Stored & operator()(const Stored & insert,
			     const ElementType & type,
			     const GhostType & ghost_type = _not_ghost);

  void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  class type_iterator : std::iterator<std::forward_iterator_tag, const ElementType> {
  public:
    typedef const ElementType   value_type;
    typedef const ElementType*  pointer;
    typedef const ElementType&  reference;
  protected:
    typedef typename ByElementType<Stored>::DataMap::const_iterator DataMapIterator;
  public:
    type_iterator(DataMapIterator & list_begin,
		  DataMapIterator & list_end,
		  UInt dim);

    type_iterator(const type_iterator & it);

    inline reference operator*();
    inline type_iterator & operator++();
    type_iterator operator++(int);
    inline bool operator==(const type_iterator & other);
    inline bool operator!=(const type_iterator & other);

  private:
    DataMapIterator list_begin;
    DataMapIterator list_end;
    UInt dim;
  };

  inline type_iterator firstType(UInt dim = 0, GhostType ghost_type = _not_ghost) const;
  inline type_iterator lastType(UInt dim = 0, GhostType ghost_type = _not_ghost) const;


protected:
  inline DataMap & getData(GhostType ghost_type);
  inline const DataMap & getData(GhostType ghost_type) const;

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

  inline Vector<T> & alloc(UInt size,
			   UInt nb_component,
			   const ElementType & type,
			   const GhostType & ghost_type);

  inline void alloc(UInt size,
		    UInt nb_component,
		    const ElementType & type);

  inline const Vector<T> & operator()(const ElementType & type,
				      const GhostType & ghost_type = _not_ghost) const;

  inline Vector<T> & operator()(const ElementType & type,
				const GhostType & ghost_type = _not_ghost);

  inline void setVector(const ElementType & type,
			const GhostType & ghost_type,
			const Vector<T> & vect);

  inline void free();
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
/* Mesh                                                                       */
/* -------------------------------------------------------------------------- */

/**
 * @class  Mesh this  contain the  coordinates of  the nodes  in  the Mesh.nodes
 * Vector,  and the  connectivity. The  connectivity  are stored  in by  element
 * types.
 *
 * To  know  all  the  element  types   present  in  a  mesh  you  can  get  the
 * Mesh::ConnectivityTypeList
 *
 * In order to loop on all element you have to loop on all types like this :
 * @code
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    UInt nb_element  = mesh.getNbElement(*it);
    UInt * conn      = mesh.getConnectivity(*it).values;

    for(UInt e = 0; e < nb_element; ++e) {
      ...
    }
  }
  @endcode
 */
class Mesh : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// constructor that create nodes coordinates array
  Mesh(UInt spatial_dimension,
       const ID & id = "mesh"
,
       const MemoryID & memory_id = 0);

  /// constructor that use an existing nodes coordinates array, by knowing its ID
  Mesh(UInt spatial_dimension,
       const ID & nodes_id,
       const ID & id = "mesh",
       const MemoryID & memory_id = 0);

  /**
   * constructor that use an existing nodes coordinates
   * array, by getting the vector of coordinates
   */
  Mesh(UInt spatial_dimension,
       Vector<Real> & nodes,
       const ID & id = "mesh",
       const MemoryID & memory_id = 0);


  virtual ~Mesh();

  /// @typedef ConnectivityTypeList list of the types present in a Mesh
  typedef std::set<ElementType> ConnectivityTypeList;

  //  typedef Vector<Real> * NormalsMap[_max_element_type];

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

  /// init a by-element-type real vector with provided ids
  template<typename T>
  void initByElementTypeVector(ByElementTypeVector<T> & v,
 			       UInt nb_component, UInt size,
			       const bool & flag_nb_node_per_elem_multiply=false) const; /// @todo: think about nicer way to do it

  /// extract coordinates of nodes from an element
  __aka_inline__ void extractNodalCoordinatesFromElement(Real * local_coords,
						 UInt * connectivity,
						 UInt n_nodes);

  /// extract coordinates of nodes from a reversed element
  __aka_inline__ void extractNodalCoordinatesFromPBCElement(Real * local_coords,
						    UInt * connectivity,
						    UInt n_nodes);

  /// convert a element to a linearized element
  __aka_inline__ UInt elementToLinearized(const Element & elem);

  /// convert a linearized element to an element
  __aka_inline__ Element linearizedToElement (UInt linearized_element);

  /// update the types offsets array for the conversions
  __aka_inline__ void updateTypesOffsets(const GhostType & ghost_type);

  /// add a Vector of connectivity for the type <type>.
  __aka_inline__ void addConnecticityType(const ElementType & type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(ID, id, const ID &);

  /// get the spatial dimension of the mesh = number of component of the coordinates
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the nodes Vector aka coordinates
  AKANTU_GET_MACRO(Nodes, *nodes, const Vector<Real> &);
  /// get the number of nodes
  AKANTU_GET_MACRO(NbNodes, nodes->getSize(), UInt);

  /// get the Vector of global ids of the nodes (only used in parallel)
  AKANTU_GET_MACRO(GlobalNodesIds, *nodes_global_ids, const Vector<UInt> &);

  /// get the global id of a node
  __aka_inline__ UInt getNodeGlobalId(UInt local_id) const;

  /// get the global number of nodes
  __aka_inline__ UInt getNbGlobalNodes() const;

  /// get the nodes type Vector
  AKANTU_GET_MACRO(NodesType, *nodes_type, const Vector<Int> &);
  __aka_inline__ Int getNodeType(UInt local_id) const;


  /// say if a node is a pure ghost node
  __aka_inline__ bool isPureGhostNode(UInt n) const;

  /// say if a node is pur local or master node
  __aka_inline__ bool isLocalOrMasterNode(UInt n) const;

  __aka_inline__ bool isLocalNode(UInt n) const;
  __aka_inline__ bool isMasterNode(UInt n) const;
  __aka_inline__ bool isSlaveNode(UInt n) const;


  AKANTU_GET_MACRO(XMin, lower_bounds[0], Real);
  AKANTU_GET_MACRO(YMin, lower_bounds[1], Real);
  AKANTU_GET_MACRO(ZMin, lower_bounds[2], Real);
  AKANTU_GET_MACRO(XMax, upper_bounds[0], Real);
  AKANTU_GET_MACRO(YMax, upper_bounds[1], Real);
  AKANTU_GET_MACRO(ZMax, upper_bounds[2], Real);

  __aka_inline__ void getLowerBounds(Real * lower) const;
  __aka_inline__ void getUpperBounds(Real * upper) const;

  __aka_inline__ void getLocalLowerBounds(Real * lower) const;
  __aka_inline__ void getLocalUpperBounds(Real * upper) const;

  /// get the number of surfaces
  AKANTU_GET_MACRO(NbSurfaces, nb_surfaces, UInt);

  /// get the connectivity Vector for a given type
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Connectivity, connectivities, UInt);

  /// @todo take out this set, if mesh can read surface id
  /// set the number of surfaces
  AKANTU_SET_MACRO(NbSurfaces, nb_surfaces, UInt);

  /// get the number of element of a type in the mesh
  __aka_inline__ UInt getNbElement(const ElementType & type, const GhostType & ghost_type = _not_ghost) const;

  // /// get the number of ghost element of a type in the mesh
  // __aka_inline__ UInt getNbGhostElement(const ElementType & type) const;

  /// get the connectivity list either for the elements or the ghost elements
  __aka_inline__ const ConnectivityTypeList & getConnectivityTypeList(const GhostType & ghost_type = _not_ghost) const;

  /// get the mesh of the internal facets
  __aka_inline__ const Mesh & getInternalFacetsMesh() const;

  /// compute the barycenter of a given element
  __aka_inline__ void getBarycenter(UInt element, const ElementType & type, Real * barycenter,
			    GhostType ghost_type = _not_ghost) const;

  /// get the surface values of facets
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(SurfaceID, surface_id, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(SurfaceID, surface_id, UInt);

  /// set the int data to the surface id vectors
  __aka_inline__ void setSurfaceIDsFromIntData(std::string & data_name);

  __aka_inline__ const Vector<UInt> & getUIntData(const ElementType & el_type,
					  const std::string & data_name,
					  const GhostType & ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Wrappers on ElementClass functions                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of nodes per element for a given element type
  static __aka_inline__ UInt getNbNodesPerElement(const ElementType & type);

  /// get the number of nodes per element for a given element type considered as
  /// a first order element
  static __aka_inline__ ElementType getP1ElementType(const ElementType & type);

  /// get spatial dimension of a type of element
  static __aka_inline__ UInt getSpatialDimension(const ElementType & type);

  /// get number of facets of a given element type
  static __aka_inline__ UInt getNbFacetsPerElement(const ElementType & type);

  /// get number of facets of a given element type
  static __aka_inline__ UInt ** getFacetLocalConnectivity(const ElementType & type);

  /// get the type of the surface element associated to a given element
  static __aka_inline__ ElementType getFacetElementType(const ElementType & type);

  /// get the pointer to the list of elements for a given type
  __aka_inline__ Vector<UInt> * getReversedElementsPBCPointer(const ElementType & type);


  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  typedef ByElementTypeUInt::type_iterator type_iterator;

  __aka_inline__ type_iterator firstType(UInt dim = 0, GhostType ghost_type = _not_ghost) const {
    return connectivities.firstType(dim, ghost_type);
  }

  __aka_inline__ type_iterator lastType(UInt dim = 0, GhostType ghost_type = _not_ghost) const {
    return connectivities.lastType(dim, ghost_type);
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

  AKANTU_GET_MACRO(NodesPointer, nodes, Vector<Real> *);

  /// get a pointer to the nodes_global_ids Vector<UInt> and create it if necessary
  __aka_inline__ Vector<UInt> * getNodesGlobalIdsPointer();

  /// get a pointer to the nodes_type Vector<Int> and create it if necessary
  __aka_inline__ Vector<Int> * getNodesTypePointer();

  /// get a pointer to the connectivity Vector for the given type and create it if necessary
  __aka_inline__ Vector<UInt> * getConnectivityPointer(const ElementType & type,
					       const GhostType & ghost_type = _not_ghost);

  /// get a pointer to the internal_facets Mesh and create it if necessary
  __aka_inline__ Mesh * getInternalFacetsMeshPointer();

  // __aka_inline__ Vector<Real> * getNormalsPointer(ElementType type) const;

  /// get a pointer to the surface_id Vector for the given type and create it if necessary
  __aka_inline__ Vector<UInt> * getSurfaceIDPointer(const ElementType & type, const GhostType & ghost_type = _not_ghost);

  /// get the UIntDataMap for a given ElementType
  __aka_inline__ UIntDataMap & getUIntDataMap(const ElementType & el_type,
				      const GhostType & ghost_type = _not_ghost);

  /// get the IntDataMap pointer (moidifyable) for a given ElementType
  __aka_inline__ Vector<UInt> * getUIntDataPointer(const ElementType & el_type,
					   const std::string & data_name,
					   const GhostType & ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// id of the mesh
  ID id;

  /// array of the nodes coordinates
  Vector<Real> * nodes;

  /// global node ids
  Vector<UInt> * nodes_global_ids;

  /// node type,  -3 pure ghost, -2  master for the  node, -1 normal node,  i in
  /// [0-N] slave node and master is proc i
  Vector<Int> * nodes_type;

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

  /// internal facets mesh
  Mesh * internal_facets_mesh;

  /// types offsets
  Vector<UInt> types_offsets;

  /// list of all existing types in the mesh
  ConnectivityTypeList ghost_type_set;
  /// ghost types offsets
  Vector<UInt> ghost_types_offsets;

  /// number of surfaces present in this mesh
  UInt nb_surfaces;

  /// surface id of the surface elements in this mesh
  ByElementTypeUInt surface_id;

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



  // /// list of elements that are reversed due to pbc
  // ByElementTypeUInt reversed_elements_pbc;
  // /// direction in which pbc are to be applied
  // UInt pbc_directions[3];

  /// list of the vectors corresponding to tags in the mesh
  ByElementTypeUIntDataMap uint_data;
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

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "mesh_inline_impl.cc"
#endif

#include "by_element_type_tmpl.hh"


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Mesh & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__


#endif /* __AKANTU_MESH_HH__ */
