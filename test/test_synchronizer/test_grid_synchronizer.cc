/**
 * @file   test_grid_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Nov 25 17:00:17 2011
 *
 * @brief  test the GridSynchronizer object
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_grid_dynamic.hh"
#include "mesh.hh"
#include "grid_synchronizer.hh"
#include "mesh_partition.hh"
#include "synchronizer_registry.hh"
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

const UInt spatial_dimension = 2;

typedef std::map<std::pair<Element, Element>, Real> pair_list;

#include "test_grid_tools.hh"

static void updatePairList(const ByElementTypeReal & barycenter,
                             const SpatialGrid<Element> & grid,
                             Real radius,
                             pair_list & neighbors,
                             neighbors_map_t<spatial_dimension>::type & neighbors_map) {
  AKANTU_DEBUG_IN();

  GhostType ghost_type = _not_ghost;

  Element e;
  e.ghost_type = ghost_type;

  // generate the pair of neighbor depending of the cell_list
  ByElementTypeReal::type_iterator it        = barycenter.firstType(0, ghost_type);
  ByElementTypeReal::type_iterator last_type = barycenter.lastType(0, ghost_type);
  for(; it != last_type; ++it) {
    // loop over quad points

    e.type = *it;
    e.element = 0;

    const Array<Real> & barycenter_vect = barycenter(*it, ghost_type);
    UInt sp = barycenter_vect.getNbComponent();

    Array<Real>::const_iterator< Vector<Real> > bary =
      barycenter_vect.begin(sp);
    Array<Real>::const_iterator< Vector<Real> > bary_end =
      barycenter_vect.end(sp);

    for(;bary != bary_end; ++bary, e.element++) {
      Point<spatial_dimension> pt1(*bary);

      SpatialGrid<Element>::CellID cell_id = grid.getCellID(*bary);
      SpatialGrid<Element>::neighbor_cells_iterator first_neigh_cell =
        grid.beginNeighborCells(cell_id);
      SpatialGrid<Element>::neighbor_cells_iterator last_neigh_cell =
        grid.endNeighborCells(cell_id);

      // loop over neighbors cells of the one containing the current element
      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
        SpatialGrid<Element>::Cell::const_iterator first_neigh_el =
          grid.beginCell(*first_neigh_cell);
        SpatialGrid<Element>::Cell::const_iterator last_neigh_el =
          grid.endCell(*first_neigh_cell);

        // loop over the quadrature point in the current cell of the cell list
        for (;first_neigh_el != last_neigh_el; ++first_neigh_el){
          const Element & elem = *first_neigh_el;

	  Array<Real>::const_iterator< Vector<Real> > neigh_it =
	    barycenter(elem.type, elem.ghost_type).begin(sp);

	  const Vector<Real> & neigh_bary = neigh_it[elem.element];

          Real distance = bary->distance(neigh_bary);
          if(distance <= radius) {
            Point<spatial_dimension> pt2(neigh_bary);
            neighbors_map[pt1].push_back(pt2);

            std::pair<Element, Element> pair = std::make_pair(e, elem);
            pair_list::iterator p = neighbors.find(pair);
            if(p != neighbors.end()) {
              AKANTU_DEBUG_ERROR("Pair already registered [" << e << " " << elem << "] -> " << p->second << " " << distance);
            } else {
              neighbors[pair] = distance;
            }
          }
        }
      }
    }
  }
  AKANTU_DEBUG_OUT();
}

class TestAccessor : public DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TestAccessor(const Mesh & mesh, const ByElementTypeReal & barycenters);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Barycenter, barycenters, Real);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  virtual UInt getNbDataForElements(const Array<Element> & elements,
                                    SynchronizationTag tag) const;
  virtual void packElementData(CommunicationBuffer & buffer,
			const Array<Element> & elements,
			SynchronizationTag tag) const;
  virtual void unpackElementData(CommunicationBuffer & buffer,
				 const Array<Element> & elements,
				 SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  const ByElementTypeReal & barycenters;
  const Mesh & mesh;
};


/* -------------------------------------------------------------------------- */
/* TestSynchronizer implementation                                            */
/* -------------------------------------------------------------------------- */
TestAccessor::TestAccessor(const Mesh & mesh,
                           const ByElementTypeReal & barycenters) : barycenters(barycenters), mesh(mesh) { }

UInt TestAccessor::getNbDataForElements(const Array<Element> & elements,
                                        __attribute__ ((unused)) SynchronizationTag tag) const {
  if(elements.getSize())
    return Mesh::getSpatialDimension(elements(0).type) * sizeof(Real) * elements.getSize();
  else
    return 0;
}

void TestAccessor::packElementData(CommunicationBuffer & buffer,
				   const Array<Element> & elements,
				   __attribute__ ((unused)) SynchronizationTag tag) const {
  UInt spatial_dimension = mesh.getSpatialDimension();
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    Vector<Real> bary(this->barycenters(element.type, element.ghost_type).storage()
                        + element.element * spatial_dimension,
                        spatial_dimension);
    buffer << bary;
  }
}

void TestAccessor::unpackElementData(CommunicationBuffer & buffer,
				     const Array<Element> & elements,
				     __attribute__ ((unused)) SynchronizationTag tag) {
  UInt spatial_dimension = mesh.getSpatialDimension();
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    Vector<Real> barycenter_loc(this->barycenters(element.type,element.ghost_type).storage()
                                  + element.element * spatial_dimension,
                                  spatial_dimension);

    Vector<Real> bary(spatial_dimension);
    buffer >> bary;
    Real tolerance = 1e-15;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if(!(std::abs(bary(i) - barycenter_loc(i)) <= tolerance))
        AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
                           << element
                           << "(barycenter[" << i << "] = " << barycenter_loc(i)
                           << " and buffer[" << i << "] = " << bary(i) << ") - tag: " << tag);
    }
  }
}

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  Real radius = 0.01;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  DistributedSynchronizer * dist = NULL;

  if(prank == 0) {
    mesh.read("triangle.msh");
    MeshPartition * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    dist = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    dist = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  mesh.computeBoundingBox();

  Real lower_bounds[spatial_dimension];
  Real upper_bounds[spatial_dimension];

  mesh.getLowerBounds(lower_bounds);
  mesh.getUpperBounds(upper_bounds);

  Vector<Real> spacing(spatial_dimension);
  Vector<Real> center(spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    spacing[i] = radius * 1.2;
    center[i] = (upper_bounds[i] + lower_bounds[i]) / 2.;
  }

  SpatialGrid<Element> grid(spatial_dimension, spacing, center);

  GhostType ghost_type = _not_ghost;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);

  ByElementTypeReal barycenters("", "", 0);
  mesh.initByElementTypeArray(barycenters, spatial_dimension, spatial_dimension);

  Element e;
  e.ghost_type = ghost_type;

  for(; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it, ghost_type);
    e.type = *it;
    Array<Real> & barycenter = barycenters(*it, ghost_type);
    barycenter.resize(nb_element);

    Array<Real>::iterator< Vector<Real> > bary_it = barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh.getBarycenter(elem, *it, bary_it->storage(), ghost_type);
      e.element = elem;
      grid.insert(e, *bary_it);
      ++bary_it;
    }
  }

  std::stringstream sstr; sstr << "mesh_" << prank << ".msh";
  mesh.write(sstr.str());

  Mesh grid_mesh(spatial_dimension, "grid_mesh", 0);
  std::stringstream sstr_grid; sstr_grid << "grid_mesh_" << prank << ".msh";
  grid.saveAsMesh(grid_mesh);
  grid_mesh.write(sstr_grid.str());



  std::cout << "Pouet 1" << std::endl;

  GridSynchronizer * grid_communicator = GridSynchronizer::createGridSynchronizer(mesh, grid);

  std::cout << "Pouet 2" << std::endl;

  ghost_type = _ghost;

  it = mesh.firstType(spatial_dimension, ghost_type);
  last_type = mesh.lastType(spatial_dimension, ghost_type);
  e.ghost_type = ghost_type;
  for(; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it, ghost_type);
    e.type = *it;
    Array<Real> & barycenter = barycenters(*it, ghost_type);
    barycenter.resize(nb_element);

    Array<Real>::iterator< Vector<Real> > bary_it = barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh.getBarycenter(elem, *it, bary_it->storage(), ghost_type);
      e.element = elem;
      grid.insert(e, *bary_it);
      ++bary_it;
    }
  }

  Mesh grid_mesh_ghost(spatial_dimension, "grid_mesh_ghost", 0);
  std::stringstream sstr_gridg; sstr_gridg << "grid_mesh_ghost_" << prank << ".msh";
  grid.saveAsMesh(grid_mesh_ghost);
  grid_mesh_ghost.write(sstr_gridg.str());

  std::cout << "Pouet 3" << std::endl;


  neighbors_map_t<spatial_dimension>::type neighbors_map;
  pair_list neighbors;

  updatePairList(barycenters, grid, radius, neighbors, neighbors_map);
  pair_list::iterator nit  = neighbors.begin();
  pair_list::iterator nend = neighbors.end();

  std::stringstream sstrp; sstrp << "pairs_" << prank;
  std::ofstream fout(sstrp.str().c_str());
  for(;nit != nend; ++nit) {
    fout << "[" << nit->first.first << "," << nit->first.second << "] -> "
         << nit->second << std::endl;
  }

  std::string file = "neighbors_ref";
  std::stringstream sstrf; sstrf << file << "_" << prank;
  file = sstrf.str();

  std::ofstream nout;
  nout.open(file);
  neighbors_map_t<spatial_dimension>::type::iterator it_n = neighbors_map.begin();
  neighbors_map_t<spatial_dimension>::type::iterator end_n = neighbors_map.end();
  for(;it_n != end_n; ++it_n) {
    std::sort(it_n->second.begin(), it_n->second.end());

    std::vector< Point<spatial_dimension> >::iterator it_v = it_n->second.begin();
    std::vector< Point<spatial_dimension> >::iterator end_v = it_n->second.end();

    nout << "####" << std::endl;
    nout << it_n->second.size() << std::endl;
    nout << it_n->first << std::endl;
    nout << "#" << std::endl;
    for(;it_v != end_v; ++it_v) {
      nout << *it_v << std::endl;
    }
  }

  fout.close();

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh, barycenters);
  SynchronizerRegistry synch_registry(test_accessor);

  synch_registry.registerSynchronizer(*grid_communicator, _gst_test);

  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);

  delete grid_communicator;
  delete dist;
  akantu::finalize();

  return EXIT_SUCCESS;
}
