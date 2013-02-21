/**
 * @file   test_mesh_partitionate_scotch_advanced.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 16 09:20:24 2012
 *
 * @brief  
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "mesh.hh"
#include "fem.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "mesh_partition_scotch.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

akantu::ElementType type = akantu::_quadrangle_4;
akantu::UInt nb_quadrature_points = 4;
akantu::UInt nb_nodes_per_element = 4;

/* -------------------------------------------------------------------------- */
void getInterfaceNodePairs(akantu::Mesh & mesh, 
			   akantu::Vector<akantu::UInt> & pairs);
void getPartition(const akantu::Mesh & mesh, 
		  const akantu::MeshPartition::EdgeLoadFunctor & fctr,
		  const akantu::Vector<akantu::UInt> & pairs,
		  akantu::MeshPartition * partition,
		  akantu::Vector<double> & parts);
void dumpParaview(akantu::Mesh & mesh, std::string name,
		  double * part_1, double * part_2);


/* -------------------------------------------------------------------------- */
class DoNotCutInterfaceFunctor : public akantu::MeshPartition::EdgeLoadFunctor {
public:
  DoNotCutInterfaceFunctor(const akantu::Mesh & mesh) : mesh(mesh) {
    
  }

  virtual inline akantu::Int operator()(const akantu::Element & el1,
					const akantu::Element & el2) const {

    const akantu::Vector<akantu::UInt> & conn_1 = this->mesh.getConnectivity(el1.type);
    const akantu::Vector<akantu::UInt> & conn_2 = this->mesh.getConnectivity(el2.type);
    std::set<akantu::UInt> nodes;

    // count number of nodes
    akantu::UInt nb_npe_1 = this->mesh.getNbNodesPerElement(el1.type);
    for (akantu::UInt i=0; i<nb_npe_1; ++i) {
      nodes.insert(conn_1(el1.element,i));
    }
    akantu::UInt nb_npe_2 = this->mesh.getNbNodesPerElement(el2.type);
    for (akantu::UInt i=0; i<nb_npe_2; ++i) {
      nodes.insert(conn_2(el2.element,i));
    }
    int max_nb_element = std::max(this->mesh.getNbElement(el1.type),
				  this->mesh.getNbElement(el2.type));
    int max_nb_npe = std::max(nb_npe_1, nb_npe_2);

    // get barycenter of elements to put different weights in vert. and horiz. dir.
    akantu::Real bc_1[2];
    mesh.getBarycenter(el1.element,el1.type,bc_1);
    akantu::Real bc_2[2];
    mesh.getBarycenter(el2.element,el2.type,bc_2);
    akantu::Real bc_diff[2];
    bc_diff[0] = std::fabs(bc_1[0] - bc_2[0]);
    bc_diff[1] = std::fabs(bc_1[1] - bc_2[1]);

    // put weights according to criterion
    int weight = 1;
    if (nodes.size() == nb_npe_1 + nb_npe_2)
      weight = max_nb_element*max_nb_npe*4;
    else if (bc_diff[0] < bc_diff[1])
      weight = 5;
    
    return weight;
  }

protected:
  const akantu::Mesh & mesh;
};


/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::debug::setDebugLevel(akantu::dblDump);
  akantu::debug::setDebugLevel(akantu::dblWarning);

  int dim = 2;
  akantu::MeshIOMSH mesh_io;

  /* ---------- check if node pairs are considered with mesh_L --------- */
  akantu::Mesh mesh_l(dim);
  mesh_io.read("squares_L.msh", mesh_l);

  // get interface node pairs
  akantu::Vector<akantu::UInt> pairs_l(0,2);
  akantu::Vector<akantu::UInt> pairs_empty(0,2);
  getInterfaceNodePairs(mesh_l,pairs_l);

  akantu::MeshPartition * partition = new akantu::MeshPartitionScotch(mesh_l, mesh_l.getSpatialDimension());

  // make normal partition -> it should cut along the interface
  akantu::Vector<double> parts(0,nb_quadrature_points);
  getPartition(mesh_l, 
	       akantu::MeshPartition::ConstEdgeLoadFunctor(),
	       pairs_empty,
	       partition,
	       parts);
  double * part = parts.storage();

  // make partition with node pairs -> it should cut perpendicular to the interface
  akantu::Vector<double> parts_adv(0,nb_quadrature_points);
  getPartition(mesh_l, 
	       akantu::MeshPartition::ConstEdgeLoadFunctor(),
	       pairs_l,
	       partition,
	       parts_adv);
  double * part_adv = parts_adv.storage();
  
  // output to visualize
#ifdef AKANTU_USE_IOHELPER
  dumpParaview(mesh_l, "test-scotch-partition-mesh-L",
	       part, part_adv);
#endif //AKANTU_USE_IOHELPER

  // check
  unsigned int nb_element = mesh_l.getNbElement(type);
  akantu::Real bb_center[2];
  mesh_l.getBarycenter(0,type,bb_center);

  // define solution for part
  unsigned int top_p_v = 0;
  unsigned int bot_p_v = 1;
  if (!(bb_center[1] > 0 && parts(0,0) == top_p_v) &&
      !(bb_center[1] < 0 && parts(0,0) == bot_p_v)) {
    top_p_v = 1;
    bot_p_v = 0;
  }
  std::cout << "top part = " << top_p_v << " | bot part = " << bot_p_v << std::endl;

  // define solution for part_adv
  unsigned int left_p_v = 0;
  unsigned int right_p_v = 1;
  if (!(bb_center[0] > 0 && parts_adv(0,0) == right_p_v) &&
      !(bb_center[0] < 0 && parts_adv(0,0) == left_p_v)) {
    left_p_v = 1;
    right_p_v = 0;
  }
  std::cout << "left part = " << left_p_v << " | right part = " << right_p_v << std::endl;

  // check
  for (akantu::UInt i=1; i<nb_element; ++i) {
    mesh_l.getBarycenter(i,type,bb_center);

    if (!(bb_center[1] > 0 && parts(i,0) == top_p_v) &&
	!(bb_center[1] < 0 && parts(i,0) == bot_p_v)) {
      std::cerr << " ** ERROR with mesh-L without node pairs" << std::endl;
      return EXIT_FAILURE;
    }
    if (!(bb_center[0] > 0 && parts_adv(i,0) == right_p_v) &&
	!(bb_center[0] < 0 && parts_adv(i,0) == left_p_v)) {
      std::cerr << " ** ERROR with mesh-L with node pairs" << std::endl;
      return EXIT_FAILURE;
    } 
  }
  
  delete partition;  


  /* ---------- check if node pairs and functor are considered with mesh_H --------- */  
  akantu::Mesh mesh_h(dim,"mesh_h",1);
  mesh_io.read("squares_H.msh", mesh_h);
  akantu::Vector<akantu::UInt> pairs_h(0,2);
  getInterfaceNodePairs(mesh_h,pairs_h);
  
  partition = new akantu::MeshPartitionScotch(mesh_h, mesh_h.getSpatialDimension());

  // make normal partition -> it should cut along the interface
  getPartition(mesh_h, 
	       akantu::MeshPartition::ConstEdgeLoadFunctor(),
	       pairs_h,
	       partition,
	       parts);
  part = parts.storage();
  
  // make partition with node pairs -> it should cut perpendicular to the interface
  DoNotCutInterfaceFunctor fctr = DoNotCutInterfaceFunctor(mesh_h);
  getPartition(mesh_h, 
	       fctr,
	       pairs_h,
	       partition,
	       parts_adv);
  part_adv = parts_adv.storage();
    
  // output to visualize
#ifdef AKANTU_USE_IOHELPER
  dumpParaview(mesh_h, "test-scotch-partition-mesh-H",
	       part, part_adv);
#endif //AKANTU_USE_IOHELPER

  // check
  nb_element = mesh_h.getNbElement(type);
  mesh_h.getBarycenter(0,type,bb_center);

  // define solution for part
  top_p_v = 0;
  bot_p_v = 1;
  if (!(bb_center[1] > 0 && parts(0,0) == top_p_v) &&
      !(bb_center[1] < 0 && parts(0,0) == bot_p_v)) {
    top_p_v = 1;
    bot_p_v = 0;
  }
  std::cout << "top part = " << top_p_v << " | bot part = " << bot_p_v << std::endl;

  // define solution for part_adv
  left_p_v = 0;
  right_p_v = 1;
  if (!(bb_center[0] > 0 && parts_adv(0,0) == right_p_v) &&
      !(bb_center[0] < 0 && parts_adv(0,0) == left_p_v)) {
    left_p_v = 1;
    right_p_v = 0;
  }
  std::cout << "left part = " << left_p_v << " | right part = " << right_p_v << std::endl;

  // check
  for (akantu::UInt i=1; i<nb_element; ++i) {
    mesh_h.getBarycenter(i,type,bb_center);

    if (!(bb_center[1] > 0 && parts(i,0) == top_p_v) &&
	!(bb_center[1] < 0 && parts(i,0) == bot_p_v)) {
      std::cerr << " ** ERROR with mesh-H with node pairs and without functor" << std::endl;
      return EXIT_FAILURE;
    }
    if (!(bb_center[0] > 0 && parts_adv(i,0) == right_p_v) &&
	!(bb_center[0] < 0 && parts_adv(i,0) == left_p_v)) {
      std::cerr << " ** ERROR with mesh-H with node pairs and with functor" << std::endl;
      return EXIT_FAILURE;
    } 
  }
  
  delete partition;  

  akantu::finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void dumpParaview(akantu::Mesh & mesh, std::string name,
		  double * part_1, double * part_2) {

  unsigned int nb_nodes = mesh.getNbNodes();
  unsigned int nb_element = mesh.getNbElement(type);

  iohelper::DumperParaview dumper;
  dumper.setMode(iohelper::TEXT);
  dumper.setPoints(mesh.getNodes().values, mesh.getSpatialDimension(), 
		   nb_nodes, name);
  dumper.setConnectivity((int*) mesh.getConnectivity(type).values,
   			 iohelper::QUAD1, nb_element, iohelper::C_MODE);

  dumper.addElemDataField("partition_1", part_1, 1, nb_element);
  dumper.addElemDataField("partition_2", part_2, 1, nb_element);
  dumper.setPrefix("paraview");
  dumper.init();
  dumper.dump();
}

/* -------------------------------------------------------------------------- */
void getPartition(const akantu::Mesh & mesh, 
		  const akantu::MeshPartition::EdgeLoadFunctor & fctr,
		  const akantu::Vector<akantu::UInt> & pairs,
		  akantu::MeshPartition * partition,
		  akantu::Vector<double> & parts) {

  //  partition->partitionate(2, akantu::MeshPartition::ConstEdgeLoadFunctor(), pairs_l);
  partition->partitionate(2, fctr, pairs);

  // get partition value for each element
  unsigned int nb_element = mesh.getNbElement(type);

  parts.resize(nb_element);
  //  double * part = new double[nb_element*nb_quadrature_points];
  akantu::UInt * part_val = partition->getPartition(type).values;
  for (unsigned int i = 0; i < nb_element; ++i)
    for (unsigned int q = 0; q < nb_quadrature_points; ++q)
      //part[i*nb_quadrature_points + q] = part_val[i];
      parts(i,q) = part_val[i];
}


/* -------------------------------------------------------------------------- */
void getInterfaceNodePairs(akantu::Mesh & mesh, 
			   akantu::Vector<akantu::UInt> & pairs) {

  // put correct number of surfaces (gmsh starts with 1 but we need 0)
  akantu::Vector<unsigned int> & surf =  const_cast<akantu::Vector<unsigned int> &>(mesh.getUIntData(akantu::_segment_2, "tag_1"));
  akantu::Vector<unsigned int>::iterator<> it = surf.begin();
  akantu::Vector<unsigned int>::iterator<> end = surf.end();
  for (;it != end; ++it)  --(*it);

  // set surface id
  mesh.setSurfaceIDsFromIntData("tag_1");
  akantu::CSR<akantu::UInt> all_surface_nodes;
  akantu::MeshUtils::buildNodesPerSurface(mesh, all_surface_nodes);

  // find top interface nodes
  akantu::Vector<akantu::UInt> top_nodes(0);
  int top = 1;
  for (akantu::CSR<akantu::UInt>::iterator node = all_surface_nodes.begin(top);
       node != all_surface_nodes.end(top);
       ++node) {
    top_nodes.push_back(*node);
  }
  
  // find bottom interface nodes
  akantu::Vector<akantu::UInt> bot_nodes(0);
  int bot = 4;
  for (akantu::CSR<akantu::UInt>::iterator node = all_surface_nodes.begin(bot);
       node != all_surface_nodes.end(bot);
       ++node) {
    bot_nodes.push_back(*node);
  }

  // verify that there is the same number of top and bottom nodes
  int nb_pairs = 0;
  if (bot_nodes.getSize() != top_nodes.getSize()) {
    std::cerr << " ** ERROR: Interface does not have the same number of top and bottom nodes" << std::endl;
  }
  else {
    nb_pairs = top_nodes.getSize();
  }

  // make pairs according to their x coordinate
  const akantu::Vector<akantu::Real> & coords = mesh.getNodes();
  int dir = 0;
  for (int i=0; i<nb_pairs; ++i) {
    int top_min = -1;
    int bot_min = -1;
    double top_min_v = std::numeric_limits<double>::max();
    double bot_min_v = std::numeric_limits<double>::max();
    for (akantu::UInt j=0; j<top_nodes.getSize(); ++j) {
      if (coords(top_nodes(j),dir) < top_min_v) {
	top_min = j;
	top_min_v = coords(top_nodes(j),dir);
      }
      if (coords(bot_nodes(j),dir) < bot_min_v) {
	bot_min = j;
	bot_min_v = coords(bot_nodes(j),dir);
      }
    }
    akantu::UInt pair[2]; 
    pair[0] = top_nodes(top_min);
    pair[1] = bot_nodes(bot_min);
    pairs.push_back(pair);
    top_nodes.erase(top_min);
    bot_nodes.erase(bot_min);
  }

}
