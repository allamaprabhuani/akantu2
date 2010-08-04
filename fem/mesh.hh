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
 */

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MESH_HH__
#define __AKANTU_MESH_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

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
       const Vector<Real> & nodes,
       const MeshID & id = "mesh",
       const MemoryID & memory_id = 0);


  virtual ~Mesh();

  typedef std::set<ElementType> ConnectivityTypeList;

  typedef Vector<UInt> * ConnectivityMap[_max_element_type];
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  AKANTU_GET_MACRO(TypeList, type_set, const ConnectivityTypeList &);

  inline Vector<Real> & getNodes() const;

  inline UInt getNbNodes() const;

  Vector<UInt> & createConnectivity(ElementType type, UInt nb_element);

  inline Vector<UInt> & getConnectivity(ElementType type) const;

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
};


/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "mesh_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Mesh & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__


#endif /* __AKANTU_MESH_HH__ */
