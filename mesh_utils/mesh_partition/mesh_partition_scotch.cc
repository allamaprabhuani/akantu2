/**
 * @file   mesh_partition_scotch.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Aug 13 11:54:11 2010
 *
 * @brief  implementation of the MeshPartitionScotch class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <cstdio>
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

#ifdef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    /// save initial graph
    FILE *fgraphinit=fopen("GraphIniFile.grf", "w");
    SCOTCH_graphSave(&scotch_graph,fgraphinit);
    fclose(fgraphinit);

    /// write geometry file
    FILE *fgeominit=fopen("GeomIniFile.xyz", "w");
    fprintf(fgeominit,"%d\n%d\n",spatial_dimension,vertnbr);

    const Vector<Real> & nodes = mesh.getNodes();

    const Mesh::ConnectivityTypeList & f_type_list = mesh.getConnectivityTypeList();
    Mesh::ConnectivityTypeList::const_iterator f_it;
    UInt out_linerized_el = 0;
    for(f_it = f_type_list.begin(); f_it != f_type_list.end(); ++f_it) {
      ElementType type = *f_it;
      if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

      UInt nb_element = mesh.getNbElement(*it);
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

	fprintf(fgeominit, "%d ", out_linerized_el++);
	for (UInt s = 0; s < spatial_dimension; ++s)
	  fprintf(fgeominit, "%f ",mid[s]);

	fprintf(fgeominit, "\n");
      }
    }
    fclose(fgeominit);
  }
#endif

  /// Partition the mesh
  SCOTCH_graphPart(&scotch_graph, nb_part, &scotch_strat, parttab);

  /// Check the graph
  AKANTU_DEBUG_ASSERT(SCOTCH_graphCheck(&scotch_graph) == 0,
		      "Partitioned graph is not consistent");

#ifdef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    /// save the partitioned graph
    FILE *fgraph=fopen("GraphFile.grf", "w");
    SCOTCH_graphSave(&scotch_graph, fgraph);
    fclose(fgraph);

    /// save the partition map
    FILE *fmap=fopen("MapFile.map", "w");
    fprintf(fmap,"%d\n",vertnbr);
    for (Int i = 0; i < vertnbr; i++) fprintf(fmap, "%d    %d\n", i, parttab[i]);
    fclose(fmap);
  }
#endif

  /// free the scotch data structures
  SCOTCH_stratExit(&scotch_strat);
  SCOTCH_graphExit(&scotch_graph);

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  UInt linerized_el = 0;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

    UInt nb_element = mesh.getNbElement(*it);

    std::stringstream sstr;    sstr    << mesh.getID() << ":partition:" << type;
    std::stringstream sstr_gi; sstr_gi << mesh.getID() << ":ghost_partition_offset:" << type;
    std::stringstream sstr_g;  sstr_g  << mesh.getID() << ":ghost_partition:" << type;

    partitions[type] = &(alloc<UInt>(sstr.str(), nb_element, 1, 0));

    ghost_partitions_offset[type] = &(alloc<UInt>(sstr_gi.str(), nb_element + 1, 1, 0));
    ghost_partitions       [type] = &(alloc<UInt>(sstr_g.str(), 0, 1, 0));

    for (UInt el = 0; el < nb_element; ++el, ++linerized_el) {
      UInt part = parttab[linerized_el];

      partitions[type]->values[el] = part;

      for (Int adj = dxadj.values[linerized_el]; adj < dxadj.values[linerized_el + 1]; ++adj) {
	UInt adj_el = dadjncy.values[adj];
	UInt adj_part = parttab[adj_el];
	if(part != adj_part) {
	  ghost_partitions[type]->push_back(adj_part);
	  ghost_partitions_offset[type]->values[el]++;
	}
      }
    }

    /// convert the ghost_partitions_offset array in an offset array
    for (UInt i = 1; i < nb_element; ++i)
      ghost_partitions_offset[type]->values[i] += ghost_partitions_offset[type]->values[i-1];
    for (UInt i = nb_element; i > 0; --i)
      ghost_partitions_offset[type]->values[i]  = ghost_partitions_offset[type]->values[i-1];
    ghost_partitions_offset[type]->values[0] = 0;
  }

  delete [] parttab;
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
