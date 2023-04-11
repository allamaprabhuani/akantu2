/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_IGFEM_ENRICHMENT_HH_
#define AKANTU_IGFEM_ENRICHMENT_HH_

#include "mesh_igfem_spherical_growing_gel.hh"
#include "mesh_sphere_intersector.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

class IGFEMEnrichment {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  IGFEMEnrichment(Mesh & mesh);
  virtual ~IGFEMEnrichment(){};

private:
  typedef std::list<Spherical::Sphere_3> Geometry;
  typedef std::map<std::string, Geometry *> GeometryMap;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// create the mesh primitives
  void initialize();

  /// get a geometry from the geometry map
  virtual Geometry & getGeometry(ID & domain) const;

  /// detect the interface
  virtual void update(ID domain = "");

  /// remove geometry
  virtual void unRegisterGeometryObject(const ID & domain = "");

  /// insert new geometry
  virtual void registerGeometryObject(Geometry & geometry,
                                      const ID & domain = "");

  /// check if a point is in a given domain
  inline bool isInside(const Vector<Real> & point, ID domain = "") const;

  /// move the interface, in this case grow the gel pockets
  virtual void moveInterface(Real new_position, ID domain = "");

  /* --------------------------------------------------------------------------
   */
  /* Accessors */
  /* --------------------------------------------------------------------------
   */
public:
  UInt getNbStandardNodes() {
    return this->intersector_sphere.getNbStandardNodes();
  }

  UInt getNbEnrichedNodes() {
    return this->intersector_sphere.getNbEnrichedNodes();
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  MeshIgfemSphericalGrowingGel<2> intersector_sphere;

  GeometryMap geometries;

  /// default geometry object
  std::string default_geometry;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "igfem_enrichment_inline_impl.hh"

} // namespace akantu
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_IGFEM_ENRICHMENT_HH_ */
