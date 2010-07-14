/**
 * @file   mesh_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 14 23:58:08 2010
 *
 * @brief  Implementation of the inline functions of the mesh class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
inline Vector<double> & Mesh::getNodes() const {
  return *nodes;
}

/* -------------------------------------------------------------------------- */
inline Vector<int> & Mesh::getConnectivity(ElementType type) const {
  AKANTU_DEBUG_IN();
  ConnectivityMap::const_iterator it = connectivities.find(type);

  AKANTU_DEBUG_ASSERT((it != connectivities.end()),
		      "The mesh " << id << " as no element of kind : "<< type);

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline Vector<int> * Mesh::getConnectivityPointer(ElementType type) const {
  AKANTU_DEBUG_IN();
  Vector<int> * conn = NULL;

  ConnectivityMap::const_iterator it = connectivities.find(type);
  if(it != connectivities.end()) {
    conn = it->second;
  }

  AKANTU_DEBUG_OUT();
  return conn;
}

/* -------------------------------------------------------------------------- */
inline const Mesh::ConnectivityMap & Mesh::getConnectivityMap() const {
  return connectivities;
}
