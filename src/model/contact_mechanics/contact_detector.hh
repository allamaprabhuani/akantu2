/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_bbox.hh"
#include "aka_common.hh"
#include "aka_grid_dynamic.hh"
#include "contact_element.hh"
#include "element_class.hh"
#include "element_group.hh"
#include "fe_engine.hh"
#include "geometry_utils.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "parsable.hh"
#include "surface_selector.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONTACT_DETECTOR_HH_
#define AKANTU_CONTACT_DETECTOR_HH_

namespace akantu {

enum class Surface { master, slave };

/* -------------------------------------------------------------------------- */

class ContactDetector : public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ContactDetector(Mesh & /*mesh*/, const ID & id = "contact_detector");

  ContactDetector(Mesh & /*mesh*/, Array<Real> positions,
                  const ID & id = "contact_detector");

  ~ContactDetector() override = default;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// performs all search steps
  void search(Array<ContactElement> & elements, Array<Real> & gaps,
              Array<Real> & normals, Array<Real> & tangents,
              Array<Real> & projections);

  /// performs global spatial search to construct spatial grids
  std::pair<BBox, std::unique_ptr<SpatialGrid<Idx>>> globalSearch() const;

  ///  performs local search to find closet master node to a slave node
  void localSearch(const BBox & intersection,
                   const SpatialGrid<Idx> & master_grid);

  /// create contact elements
  void createContactElements(Array<ContactElement> & elements,
                             Array<Real> & gaps, Array<Real> & normals,
                             Array<Real> & tangents, Array<Real> & projections);

private:
  /// reads the input file to get contact detection options
  void parseSection(const ParserSection & section) override;

  /* ------------------------------------------------------------------------ */
  /* Inline Methods                                                           */
  /* ------------------------------------------------------------------------ */
public:
  /// checks whether the natural projection is valid or not
  template <class Derived,
            std::enable_if_t<aka::is_vector_v<Derived>> * = nullptr>
  inline bool
  checkValidityOfProjection(Eigen::MatrixBase<Derived> & projection) const;

  /// extracts the coordinates of an element
  template <class Derived>
  inline void coordinatesOfElement(const Element & el,
                                   Eigen::MatrixBase<Derived> & coords) const;

  /// computes the optimal cell size for grid
  template <class Derived,
            std::enable_if_t<aka::is_vector_v<Derived>> * = nullptr>
  inline void computeCellSpacing(Eigen::MatrixBase<Derived> & spacing) const;

  /// constructs a grid containing nodes lying within bounding box
  inline void constructGrid(SpatialGrid<Idx> & grid, BBox & bbox,
                            const Array<Idx> & nodes_list) const;

  /// constructs the bounding box based on nodes list
  inline void constructBoundingBox(BBox & bbox,
                                   const Array<Idx> & nodes_list) const;

  /// computes the maximum in radius for a given mesh
  inline void computeMaximalDetectionDistance();

  /// constructs the connectivity for a contact element
  inline Vector<Idx> constructConnectivity(Idx & slave,
                                           const Element & master) const;

  /// computes normal on an element
  template <class Derived,
            std::enable_if_t<aka::is_vector_v<Derived>> * = nullptr>
  inline void computeNormalOnElement(const Element & element,
                                     Eigen::MatrixBase<Derived> & normal) const;

  /// extracts vectors which forms the plane of element
  template <class Derived>
  inline void vectorsAlongElement(const Element & el,
                                  Eigen::MatrixBase<Derived> & vectors) const;

  /// computes the gap between slave and its projection on master
  /// surface
  template <
      class Derived1, class Derived2,
      std::enable_if_t<aka::are_vectors<Derived1, Derived2>::value> * = nullptr>
  inline Real computeGap(const Eigen::MatrixBase<Derived1> & slave,
                         const Eigen::MatrixBase<Derived2> & master) const;

  /// filter boundary elements
  inline void filterBoundaryElements(const Array<Element> & elements,
                                     Array<Element> & boundary_elements) const;

  /// checks whether self contact condition leads to a master element
  /// which is closet but not orthogonally opposite to slave surface
  template <class Derived,
            std::enable_if_t<aka::is_vector_v<Derived>> * = nullptr>
  inline bool
  isValidSelfContact(const Idx & slave_node, const Real & gap,
                     const Eigen::MatrixBase<Derived> & normal) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the mesh
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &)

  /// returns the maximum detection distance
  AKANTU_GET_MACRO(MaximumDetectionDistance, max_dd, Real);
  AKANTU_SET_MACRO(MaximumDetectionDistance, max_dd, Real);

  /// returns the bounding box extension
  AKANTU_GET_MACRO(MaximumBoundingBox, max_bb, Real);
  AKANTU_SET_MACRO(MaximumBoundingBox, max_bb, Real);

  /// returns the minimum detection distance
  AKANTU_GET_MACRO(MinimumDetectionDistance, min_dd, Real);
  AKANTU_SET_MACRO(MinimumDetectionDistance, min_dd, Real);

  AKANTU_GET_MACRO_NOT_CONST(Positions, positions, Array<Real> &);
  AKANTU_SET_MACRO(Positions, positions, Array<Real>);

  AKANTU_GET_MACRO_NOT_CONST(SurfaceSelector, *surface_selector,
                             SurfaceSelector &);
  AKANTU_SET_MACRO(SurfaceSelector, surface_selector,
                   std::shared_ptr<SurfaceSelector>);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// maximal detection distance for grid spacing
  Real max_dd;

  /// minimal detection distance for grid spacing
  Real min_dd;

  /// maximal bounding box extension
  Real max_bb;

  /// tolerance for finding natural projection
  Real projection_tolerance;

  /// iterations for finding natural projection
  Int max_iterations;

  /// tolerance for extending a master elements on all sides
  Real extension_tolerance;

  /// Mesh
  Mesh & mesh;

  /// dimension of the model
  Int spatial_dimension{0};

  /// node selector for selecting master and slave nodes
  std::shared_ptr<SurfaceSelector> surface_selector;

  /// contact pair slave node to closet master node
  std::vector<std::pair<Idx, Idx>> contact_pairs;

  /// contains the updated positions of the nodes
  Array<Real> positions;

  /// type of detection explicit/implicit
  DetectionType detection_type;
};

} // namespace akantu

#include "contact_detector_inline_impl.hh"

#endif /* AKANTU_CONTACT_DETECTOR_HH_ */
