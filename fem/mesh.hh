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
/// to store data Vector<Real> by element type
typedef Vector<Real> * ByElementTypeReal[_max_element_type];

/// to store data Vector<Int> by element type
typedef Vector<Int>  * ByElementTypeInt [_max_element_type];

/// to store data Vector<UInt> by element type
typedef Vector<UInt> * ByElementTypeUInt[_max_element_type];

/* -------------------------------------------------------------------------- */
class Element {
public:
  Element(ElementType type = _not_defined, UInt element = 0) :
    type(type), element(element) {};

  Element(const Element & element) {
    this->type    = element.type;
    this->element = element.element;
  }

  ~Element() {};

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
public:
  ElementType type;
  UInt element;
};

/* -------------------------------------------------------------------------- */

/**
 * @class Mesh this contain the coordinates of the nodes in the Mesh.nodes Vector,
 * and the connectivity. The connectivity are stored in by element types.
 *
 * To know all the element types present in a mesh you can get the Mesh::ConnectivityTypeList
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
       const MeshID & id = "mesh",
       const MemoryID & memory_id = 0);

  /// constructor that use an existing nodes coordinates array, by knowing its ID
  Mesh(UInt spatial_dimension,
       const VectorID & nodes_id,
       const MeshID & id = "mesh",
       const MemoryID & memory_id = 0);

  /**
   * constructor that use an existing nodes coordinates
   * array, by getting the vector of coordinates
   */
  Mesh(UInt spatial_dimension,
       Vector<Real> & nodes,
       const MeshID & id = "mesh",
       const MemoryID & memory_id = 0);


  virtual ~Mesh();

  /// @typedef ConnectivityTypeList list of the types present in a Mesh
  typedef std::set<ElementType> ConnectivityTypeList;

  typedef Vector<Real> * NormalsMap[_max_element_type];

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
  void initByElementTypeRealVector(ByElementTypeReal & v, UInt nb_component,
				   UInt dimension,
				   const std::string & obj_id,
				   const std::string & vec_id,
				   GhostType ghost_type);

  /// init a by-element-type real vector with provided ids
  void initByElementTypeUIntVector(ByElementTypeUInt & v,UInt nb_component,
				   UInt dimension,
				   const std::string & obj_id,
				   const std::string & vec_id,
				   GhostType ghost_type=_not_ghost);

  /// extract coordinates of nodes from an element
  inline void extractNodalCoordinatesFromElement(Real * local_coords,
						 UInt * connectivity,
						 UInt n_nodes);

  /// convert a element to a linearized element
  inline UInt elementToLinearized(const Element & elem);

  /// convert a linearized element to an element
  inline Element linearizedToElement (UInt linearized_element);

  /// update the types offsets array for the conversions
  inline void updateTypesOffsets();

  /// convert a element to a linearized element
  inline UInt ghostElementToLinearized(const Element & elem);

  /// convert a linearized element to an element
  inline Element ghostLinearizedToElement (UInt linearized_element);

  /// update the types offsets array for the conversions
  inline void updateGhostTypesOffsets();

  /// add a Vector of connectivity for the type <type>.
  inline void addConnecticityType(const ElementType & type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(ID, id, const MeshID &);

  /// get the spatial dimension of the mesh = number of component of the coordinates
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the nodes Vector aka coordinates
  AKANTU_GET_MACRO(Nodes, *nodes, const Vector<Real> &);
  /// get the number of nodes
  AKANTU_GET_MACRO(NbNodes, nodes->getSize(), UInt);

  /// get the Vector of global ids of the nodes (only used in parallel)
  AKANTU_GET_MACRO(GlobalNodesIds, *nodes_global_ids, const Vector<UInt> &);

  /// get the global id of a node
  inline UInt getNodeGlobalId(UInt local_id) const;

  /// get the global number of nodes
  inline UInt getNbGlobalNodes() const;

  /// get the nodes type Vector
  AKANTU_GET_MACRO(NodesType, *nodes_type, const Vector<Int> &);

  /// say if a node is local or not
  inline bool isLocalNode(UInt n) const;

  /// say if a node is pur local or master node
  inline bool isLocalOrMasterNode(UInt n) const;

  AKANTU_GET_MACRO(XMin, xmin[0], UInt);
  AKANTU_GET_MACRO(YMin, xmin[1], UInt);
  AKANTU_GET_MACRO(ZMin, xmin[2], UInt);

  AKANTU_GET_MACRO(XMax, xmax[0], UInt);
  AKANTU_GET_MACRO(YMax, xmax[1], UInt);
  AKANTU_GET_MACRO(ZMax, xmax[2], UInt);

  /// get the number of surfaces
  AKANTU_GET_MACRO(NbSurfaces, nb_surfaces, UInt);

  /// get the connectivity Vector for a given type
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Connectivity, connectivities, const Vector<UInt> &);
  /// get the connecticity of ghost elements of a given type
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(GhostConnectivity, ghost_connectivities, const Vector<UInt> &);

  //  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Normals, normals, const Vector<Real> &);

  /// @todo take out this set, if mesh can read surface id
  /// set the number of surfaces
  AKANTU_SET_MACRO(NbSurfaces, nb_surfaces, UInt);

  /// get the number of element of a type in the mesh
  inline UInt getNbElement(const ElementType & type) const;

  /// get the number of ghost element of a type in the mesh
  inline UInt getNbGhostElement(const ElementType & type) const;

  /// get the connectivity list either for the elements or the ghost elements
  inline const ConnectivityTypeList & getConnectivityTypeList(GhostType ghost_type = _not_ghost) const;

  /// get the mesh of the internal facets
  inline const Mesh & getInternalFacetsMesh() const;

  /// compute the barycenter of a given element
  inline void getBarycenter(UInt element, ElementType type, Real * barycenter,
			    GhostType ghost_type = _not_ghost) const;

  /// get the surface values of facets
  inline const Vector<UInt> & getSurfaceId(const ElementType & type) const;

  /// set the int data to the surface id vectors
  inline void setSurfaceIdsFromIntData(std::string & data_name);

  inline const Vector<UInt> & getUIntData(const ElementType & el_type,
					  const std::string & data_name) const;

  /* ------------------------------------------------------------------------ */
  /* Wrappers on ElementClass functions                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of nodes per element for a given element type
  static inline UInt getNbNodesPerElement(const ElementType & type);

  /// get the number of nodes per element for a given element type considered as
  /// a first order element
  static inline ElementType getP1ElementType(const ElementType & type);

  /// get spatial dimension of a type of element
  static inline UInt getSpatialDimension(const ElementType & type);

  /// get number of facets of a given element type
  static inline UInt getNbFacetsPerElement(const ElementType & type);

  /// get number of facets of a given element type
  static inline UInt ** getFacetLocalConnectivity(const ElementType & type);

  /// get the type of the surface element associated to a given element
  static inline ElementType getFacetElementType(const ElementType & type);

private:
  friend class MeshIOMSH;
  friend class MeshIODiana;
  friend class MeshUtils;
  friend class Communicator;

  AKANTU_GET_MACRO(NodesPointer, nodes, Vector<Real> *);

  /// get a pointer to the nodes_global_ids Vector<UInt> and create it if necessary
  inline Vector<UInt> * getNodesGlobalIdsPointer();

  /// get a pointer to the nodes_type Vector<Int> and create it if necessary
  inline Vector<Int> * getNodesTypePointer();

  /// get a pointer to the connectivity Vector for the given type and create it if necessary
  inline Vector<UInt> * getConnectivityPointer(ElementType type);

  /// get a pointer to the ghost_connectivity Vector for the given type and create it if necessary
  inline Vector<UInt> * getGhostConnectivityPointer(ElementType type);

  /// get a pointer to the internal_facets Mesh and create it if necessary
  inline Mesh * getInternalFacetsMeshPointer();

  // inline Vector<Real> * getNormalsPointer(ElementType type) const;

  /// get a pointer to the surface_id Vector for the given type and create it if necessary
  inline Vector<UInt> * getSurfaceIdPointer(ElementType type);

  typedef std::map<std::string, Vector<UInt> * > UIntDataMap;
  typedef UIntDataMap ByElementTypeUIntDataMap[_max_element_type];
  /// get the IntDataMap for a given ElementType
  inline UIntDataMap & getUIntDataMap(const ElementType & el_type) { return uint_data[el_type]; };
  inline UIntDataMap & getGhostUIntDataMap(const ElementType & el_type) { return ghost_uint_data[el_type]; };

  inline Vector<UInt> * getUIntDataPointer(const ElementType & el_type, const std::string & data_name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// id of the mesh
  MeshID id;

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

  /// all class of elements present in this mesh (for heterogenous meshes)
  ByElementTypeUInt ghost_connectivities;

  /// list of all existing types in the mesh
  ConnectivityTypeList ghost_type_set;

  /// ghost types offsets
  Vector<UInt> ghost_types_offsets;

  /// number of surfaces present in this mesh
  UInt nb_surfaces;

  /// surface id of the surface elements in this mesh
  ByElementTypeUInt surface_id;
  ByElementTypeUInt ghost_surface_id;

  /// min of coordinates
  Real xmin[3];
  /// max of coordinates
  Real xmax[3];

  /// list of elements that are reversed due to pbc
  ByElementTypeUInt reversed_elements_pbc;

  /// list of the vectors corresponding to tags in the mesh
  ByElementTypeUIntDataMap uint_data;
  ByElementTypeUIntDataMap ghost_uint_data;
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

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Mesh & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__


#endif /* __AKANTU_MESH_HH__ */
