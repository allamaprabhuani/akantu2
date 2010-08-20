/**
 * @file   mesh.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 11:53:53 2010
 *
 * @brief  the class representing the meshes
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 * @todo struct element + function linerized index <-> element
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

class Element {
public:
  Element(ElementType type = _not_defined, UInt element = 0) :
    type(type), element(element) {};

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
public:
  ElementType type;
  UInt element;
};

/* -------------------------------------------------------------------------- */
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

  typedef std::set<ElementType> ConnectivityTypeList;

  typedef Vector<UInt> * ConnectivityMap[_max_element_type];

private:
  /// initialize the connectivity to NULL
  void initConnectivities();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  Vector<UInt> & createConnectivity(ElementType type, UInt nb_element);

#ifdef AKANTU_USE_MPI
  Vector<UInt> & createGhostConnectivity(ElementType type, UInt nb_element);
#endif //AKANTU_USE_MPI

  /// convert a element to a linearized element
  inline UInt elementToLinearized(const Element & elem);

  /// convert a linearized element to an element
  inline Element linearizedToElement (UInt linearized_element);

  /// update the types offsets array for the conversions
  inline void updateTypesOffsets();

#ifdef AKANTU_USE_MPI
  /// convert a element to a linearized element
  inline UInt ghostElementToLinearized(const Element & elem);

  /// convert a linearized element to an element
  inline Element ghostLinearizedToElement (UInt linearized_element);

  /// update the types offsets array for the conversions
  inline void updateGhostTypesOffsets();
#endif //AKANTU_USE_MPI


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(ID, id, const MeshID &);

  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Connectivity, connectivities, Vector<UInt> &);

#ifdef AKANTU_USE_MPI
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Connectivity, ghost_connectivities, Vector<UInt> &);
#endif //AKANTU_USE_MPI

  /// get the number of element of a type in the mesh
  inline UInt getNbElement(const ElementType & type) const;

  inline const ConnectivityTypeList & getConnectivityTypeList(bool local = true) const;

  AKANTU_GET_MACRO(Nodes, *nodes, Vector<Real> &);
  AKANTU_GET_MACRO(NbNodes, nodes->getSize(), UInt);

  /// get the number of nodes per element for a given element type
  static inline UInt getNbNodesPerElement(const ElementType & type);

  /// get the number of nodes per element for a given element type considered as
  /// a first order element
  static inline UInt getNbNodesPerElementP1(const ElementType & type);

  /// get spatial dimension of a type of element
  static inline UInt getSpatialDimension(const ElementType & type);

  /// get the type of the surface element associated to a given element
  static inline const ElementType getSurfaceElementType(const ElementType & type);

private:
  friend class MeshIOMSH;

  inline Vector<UInt> * getConnectivityPointer(ElementType type) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// id of the mesh
  MeshID id;

  /// array of the nodes coordinates
  Vector<Real> * nodes;

  /// boolean to know if the nodes have to be deleted with the mesh or not
  bool created_nodes;

  /// all class of elements present in this mesh (for heterogenous meshes)
  ConnectivityMap connectivities;

  /// list of all existing types in the mesh
  ConnectivityTypeList type_set;

  /// the spatial dimension of this mesh
  UInt spatial_dimension;

  /// types offsets
  Vector<UInt> types_offsets;

#ifdef AKANTU_USE_MPI
  /// all class of elements present in this mesh (for heterogenous meshes)
  ConnectivityMap ghost_connectivities;

  /// list of all existing types in the mesh
  ConnectivityTypeList ghost_type_set;

  /// ghost types offsets
  Vector<UInt> ghost_types_offsets;
#endif //AKANTU_USE_MPI
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
