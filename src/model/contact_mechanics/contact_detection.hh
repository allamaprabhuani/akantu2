/**
 * @file contact_detection.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Sep 12 2018
 * @date last modification: Tue Oct 23 2018
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
#include "fe_engine.hh"
#include "parsable.hh"
#include "element_group.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_DETECTION_HH__
#define __AKANTU_CONTACT_DETECTION_HH__


namespace akantu { 
  
/* -------------------------------------------------------------------------- */ 

class ContactDetection :
    private Memory, public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  ContactDetection(Mesh &, std::string , std::string , const ID & id = "contact_detection",
		   UInt memory_id = 0);

  ~ContactDetection() = default;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ///
  void search();

  /// performs global spatial search
  void globalSearch();

  ///  performs local search to create contact element
  void localSearch();

  /// constructs a grid containing nodes lying within bounding box
  void constructGrid(SpatialGrid<UInt> &, BBox &, const Array<UInt> &);

  /// constructs the bounding box based on nodes list
  void constructBoundingBox(BBox &, const Array<UInt> &);

  /// computes the optimal cell size for grid 
  void computeCellSpacing(Vector<Real> &);

  /// computes the maximum in radius for a given mesh 
  void getMaximalDetectionDistance();
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// maximal detection distance
  Real max_dd;
  
  /// dimension of the model
  UInt spatial_dimension{0};

  /// Mesh
  Mesh & mesh;

  /// id for master surface/curve
  std::string master_id;

  /// id for slave surface/curve
  std::string slave_id;

};
  
} // namespace akantu


#endif /* __AKANTU_CONTACT_DETECTION_HH__ */
