/**
 * @file   mesh_partition_scotch.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Aug 13 11:54:11 2010
 *
 * @brief  implementation of the MeshPartitionScotch class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include <cstdio>
#include <fstream>
extern "C" {
#include <scotch.h>
}
/* -------------------------------------------------------------------------- */
#include "mesh_partition_scotch.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshPartitionScotch::MeshPartitionScotch(const Mesh & mesh, UInt spatial_dimension,
					 const MemoryID & memory_id) :
  MeshPartition(mesh, spatial_dimension, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartitionScotch::partitionate(UInt nb_part) {
  AKANTU_DEBUG_IN();

  nb_partitions = nb_part;

  AKANTU_DEBUG_INFO("Partitioning the mesh " << mesh.getID()
		    << " in " << nb_part << " parts.");

  Vector<Int> dxadj;
  Vector<Int> dadjncy;

  buildDualGraph(dxadj, dadjncy);

  /// variables that will hold our structures in scotch format
  SCOTCH_Graph   scotch_graph;
  SCOTCH_Strat   scotch_strat;

  /// description number and arrays for struct mesh for scotch
  SCOTCH_Num baseval = 0;                       //base numbering for element and
						//nodes (0 -> C , 1 -> fortran)
  SCOTCH_Num vertnbr = dxadj.getSize() - 1;     //number of vertexes
  SCOTCH_Num *parttab;                          //array of partitions
  SCOTCH_Num edgenbr = dxadj.values[vertnbr];   //twice  the number  of "edges"
						//(an "edge" bounds two nodes)
  SCOTCH_Num * verttab = dxadj.values;          //array of start indices in edgetab
  SCOTCH_Num * vendtab = NULL;                  //array of after-last indices in edgetab
  SCOTCH_Num * velotab = NULL;                  //integer  load  associated with
						//every vertex ( optional )
  SCOTCH_Num *edlotab = NULL;                   //integer  load  associated with
						//every edge ( optional )
  SCOTCH_Num *edgetab = dadjncy.values;         // adjacency array of every vertex
  SCOTCH_Num *vlbltab = NULL;                   // vertex label array (optional)

  /// Allocate space for Scotch arrays
  parttab = new SCOTCH_Num[vertnbr];

  /// Initialize the strategy structure
  SCOTCH_stratInit (&scotch_strat);

  /// Initialize the graph structure
  SCOTCH_graphInit(&scotch_graph);

  /// Build the graph from the adjacency arrays
  SCOTCH_graphBuild(&scotch_graph, baseval,
		    vertnbr, verttab, vendtab,
		    velotab, vlbltab, edgenbr,
		    edgetab, edlotab);

  /// Check the graph
  AKANTU_DEBUG_ASSERT(SCOTCH_graphCheck(&scotch_graph) == 0,
		      "Graph to partition is not consistent");

#ifndef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    /// save initial graph
    FILE *fgraphinit = fopen("GraphIniFile.grf", "w");
    SCOTCH_graphSave(&scotch_graph,fgraphinit);
    fclose(fgraphinit);

    /// write geometry file
    std::ofstream fgeominit;
    fgeominit.open("GeomIniFile.xyz");
    fgeominit << spatial_dimension << std::endl << vertnbr << std::endl;

    const Vector<Real> & nodes = mesh.getNodes();

    const Mesh::ConnectivityTypeList & f_type_list = mesh.getConnectivityTypeList();
    Mesh::ConnectivityTypeList::const_iterator f_it;
    UInt out_linerized_el = 0;
    for(f_it = f_type_list.begin(); f_it != f_type_list.end(); ++f_it) {
      ElementType type = *f_it;
      if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

      UInt nb_element = mesh.getNbElement(*f_it);
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
      const Vector<UInt> & connectivity = mesh.getConnectivity(type);

      Real mid[spatial_dimension] ;
      for (UInt el = 0; el < nb_element; ++el) {
	memset(mid, 0, spatial_dimension*sizeof(Real));
	for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	  UInt node = connectivity.values[nb_nodes_per_element * el + n];
	  for (UInt s = 0; s < spatial_dimension; ++s)
	    mid[s] += ((Real) nodes.values[node * spatial_dimension + s]) / ((Real) nb_nodes_per_element);
	}

	fgeominit << out_linerized_el++ << " ";
	for (UInt s = 0; s < spatial_dimension; ++s)
	  fgeominit << mid[s] << " ";

	fgeominit << std::endl;;
      }
    }
    fgeominit.close();
  }
#endif

  /// Partition the mesh
  SCOTCH_graphPart(&scotch_graph, nb_part, &scotch_strat, parttab);

  /// Check the graph
  AKANTU_DEBUG_ASSERT(SCOTCH_graphCheck(&scotch_graph) == 0,
		      "Partitioned graph is not consistent");

#ifndef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    /// save the partitioned graph
    FILE *fgraph=fopen("GraphFile.grf", "w");
    SCOTCH_graphSave(&scotch_graph, fgraph);
    fclose(fgraph);

    /// save the partition map
    std::ofstream fmap;
    fmap.open("MapFile.map");
    fmap << vertnbr << std::endl;
    for (Int i = 0; i < vertnbr; i++) fmap << i << "    " << parttab[i] << std::endl;
    fmap.close();
  }
#endif

  /// free the scotch data structures
  SCOTCH_stratExit(&scotch_strat);
  SCOTCH_graphFree(&scotch_graph);
  SCOTCH_graphExit(&scotch_graph);

  fillPartitionInformations(mesh, parttab);

  delete [] parttab;
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
