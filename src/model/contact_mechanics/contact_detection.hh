/**
 * @file contact_detection.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Sep 12 2018
 * @date last modification: Fri Sep 21 2018
 *
 * @brief  Mother class for all detection algorithms
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "aka_grid_dynamic.hh"
#include "aka_bbox.hh"
#include "mesh.hh"
#include "mesh_io.hh"
/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_CONTACT_DETECTION_HH__
#define __AKANTU_CONTACT_DETECTION_HH__


namespace akantu {

struct NodeInfo {
    NodeInfo() {}
    NodeInfo(UInt spatial_dimension) : position(spatial_dimension) {}
    NodeInfo(UInt node, const Vector<Real> & position)
        : node(node), position(position) {
    }

    NodeInfo(const NodeInfo & other)
        : node(other.node), position(other.position)
    {}

    UInt node{0};
    Vector<Real> position;
  };
  
  

  
class ContactDetection {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ContactDetection();

  ContactDetection(Mesh &);

  ~ContactDetection() = default;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ///
  virtual void search();

private:
  ///
  void globalSearch();

  ///
  void localSearch();

  ///
  void constructGrid(SpatialGrid<Element> &);

  ///
  void constructBoundingBox(BBox &, const Array<UInt> &);

  ///
  void computeCellSpacing(Vector<Real> &);

  ///
  void computeMaximalDetectionDistance();
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// maximal detection distance
  Real d_max;
  
  ///
  UInt spatial_dimension{0};

  /// Mesh
  Mesh & mesh;


};
  
} // namespace akantu


#endif /* __AKANTU_CONTACT_DETECTION_HH__ */
