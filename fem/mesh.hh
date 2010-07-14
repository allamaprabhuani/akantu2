/**
 * @file   mesh.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 11:53:53 2010
 *
 * @brief  the class representing the meshes
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MESH_HH__
#define __AKANTU_MESH_HH__

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "memory.hh"
#include "vector.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class Mesh : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// constructor that create nodes coordinates array
  Mesh(unsigned int spatial_dimension,
       const MeshID & id = "mesh",
       const MemoryID & memory_id = 0);

  /// constructor that use an existing nodes coordinates array, by knowing its ID
  Mesh(unsigned int spatial_dimension,
       const VectorID & nodes_id,
       const MeshID & id = "mesh",
       const MemoryID & memory_id = 0);

  /**
   * constructor that use an existing nodes coordinates
   * array, by getting the vector of coordinates
   */
  Mesh(unsigned int spatial_dimension,
       const Vector<double> & nodes,
       const MeshID & id = "mesh",
       const MemoryID & memory_id = 0);


  virtual ~Mesh() {};

  typedef std::map<ElementType, Vector<int> *> ConnectivityMap;

  typedef std::map<ElementType, Vector<int> *>::iterator ConnectivityIterator;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline Vector<double> & getNodes() const;

  Vector<int> & createConnectivity(ElementType type, unsigned int nb_element);

  inline Vector<int> & getConnectivity(ElementType type) const;

  inline const ConnectivityMap & getConnectivityMap() const;

private:
  friend class MeshIOMSH;

  inline Vector<int> * getConnectivityPointer(ElementType type) const;


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// id of the mesh
  MeshID id;

  /// array of the nodes coordinates
  Vector<double> * nodes;

  /// all class of elements present in this mesh (for heterogenous meshes)
  ConnectivityMap connectivities;

  /// the spatial dimension of this mesh
  unsigned int spatial_dimension;
};


/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "mesh_inline_impl.cc"

__END_AKANTU__


#endif /* __AKANTU_MESH_HH__ */
