/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AKANTU_NEIGHBORHOOD_BASE_HH_
#define AKANTU_NEIGHBORHOOD_BASE_HH_
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "data_accessor.hh"
#include "grid_synchronizer.hh"
#include "integration_point.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
class Model;
template <class T> class SpatialGrid;
class RemovedElementsEvent;
} // namespace akantu

namespace akantu {
class NeighborhoodBase : public DataAccessor<Element>,
                         public SynchronizerRegistry {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NeighborhoodBase(Model & model,
                   const ElementTypeMapArray<Real> & quad_coordinates,
                   const ID & id = "neighborhood");

  using PairList = std::vector<std::pair<IntegrationPoint, IntegrationPoint>>;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  /// intialize the neighborhood
  virtual void initNeighborhood();

  /// initialize the material computed parameter
  inline void insertIntegrationPoint(const IntegrationPoint & quad,
                                     const VectorProxy<const Real> & coords);

  /// create the pairs of quadrature points
  void updatePairList();

  /// save the pairs of quadrature points in a file
  void savePairs(const std::string & filename) const;

  /// save the coordinates of all neighbors of a quad
  void saveNeighborCoords(const std::string & filename) const;

  /// create grid synchronizer and exchange ghost cells
  virtual void createGridSynchronizer() = 0;
  virtual void synchronize(DataAccessor<Element> & data_accessor,
                           const SynchronizationTag & tag) = 0;

  /// inherited function from MeshEventHandler
  virtual void onElementsRemoved(const Array<Element> & element_list,
                                 const ElementTypeMapArray<Idx> & new_numbering,
                                 const RemovedElementsEvent & event);

protected:
  /// create the grid
  void createGrid();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, Int);
  AKANTU_GET_MACRO(Model, model, const Model &);
  /// return the object handling synchronizers
  const PairList & getPairLists(GhostType type) { return pair_list[type]; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  ID id;

  /// the model to which the neighborhood belongs
  Model & model;

  /// Radius of impact: to determine if two quadrature points influence each
  /// other
  Real neighborhood_radius{0.};

  /**
   * the pairs of quadrature points
   * 0: not ghost to not ghost
   * 1: not ghost to ghost
   */
  std::map<GhostType, PairList> pair_list;

  /// the regular grid to construct/update the pair lists
  std::unique_ptr<SpatialGrid<IntegrationPoint>> spatial_grid;

  bool is_creating_grid{false};

  /// the grid synchronizer for parallel computations
  std::unique_ptr<GridSynchronizer> grid_synchronizer;

  /// the quadrature point positions
  const ElementTypeMapArray<Real> & quad_coordinates;

  /// the spatial dimension of the problem
  const Int spatial_dimension;
};

} // namespace akantu

#include "neighborhood_base_inline_impl.hh"

#endif /* AKANTU_NEIGHBORHOOD_BASE_HH_ */
