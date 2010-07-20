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
inline Vector<Real> & Mesh::getNodes() const {
  return *nodes;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodes() const {
  return nodes->getSize();
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> & Mesh::getConnectivity(ElementType type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(connectivities[type] != NULL,
		      "The mesh " << id << " as no element of kind : "<< type);

  AKANTU_DEBUG_OUT();
  return *connectivities[type];
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> * Mesh::getConnectivityPointer(ElementType type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return connectivities[type];
}

