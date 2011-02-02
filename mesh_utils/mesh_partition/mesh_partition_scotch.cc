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

/* -------------------------------------------------------------------------- */
#include "mesh_partition_scotch.hh"
#include "mesh_utils.hh"

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
SCOTCH_Mesh * MeshPartitionScotch::createMesh() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes = mesh.getNbNodes();

  UInt total_nb_element = 0;
  UInt nb_edge = 0;

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != spatial_dimension) continue;

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    total_nb_element += nb_element;
    nb_edge += nb_element * nb_nodes_per_element;
  }

  SCOTCH_Num vnodbas = 0;
  SCOTCH_Num vnodnbr = nb_nodes;

  SCOTCH_Num velmbas = vnodnbr;
  SCOTCH_Num velmnbr = total_nb_element;

  SCOTCH_Num * verttab = new SCOTCH_Num[vnodnbr + velmnbr + 1];
  SCOTCH_Num * vendtab = verttab + 1;

  SCOTCH_Num * velotab = NULL;
  SCOTCH_Num * vnlotab = NULL;
  SCOTCH_Num * vlbltab = NULL;

  memset(verttab, 0, (vnodnbr + velmnbr + 1) * sizeof(SCOTCH_Num));

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != spatial_dimension) continue;

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    const Vector<UInt> & connectivity = mesh.getConnectivity(type);

    /// count number of occurrence of each node
    for (UInt el = 0; el < nb_element; ++el) {
      UInt * conn_val = connectivity.values + el * nb_nodes_per_element;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	verttab[*(conn_val++)]++;
      }
    }
  }

  /// convert the occurrence array in a csr one
  for (UInt i = 1; i < nb_nodes; ++i) verttab[i] += verttab[i-1];
  for (UInt i = nb_nodes; i > 0; --i) verttab[i]  = verttab[i-1];
  verttab[0] = 0;

  /// rearrange element to get the node-element list
  SCOTCH_Num   edgenbr = verttab[vnodnbr] + nb_edge;
  SCOTCH_Num * edgetab = new SCOTCH_Num[edgenbr];

  UInt linearized_el = 0;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != spatial_dimension) continue;

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    const Vector<UInt> & connectivity = mesh.getConnectivity(type);

    for (UInt el = 0; el < nb_element; ++el, ++linearized_el) {
      UInt * conn_val = connectivity.values + el * nb_nodes_per_element;
      for (UInt n = 0; n < nb_nodes_per_element; ++n)
	edgetab[verttab[*(conn_val++)]++] = linearized_el + velmbas;
    }
  }

  for (UInt i = nb_nodes; i > 0; --i) verttab[i]  = verttab[i-1];
  verttab[0] = 0;

  SCOTCH_Num * verttab_tmp = verttab + vnodnbr + 1;
  SCOTCH_Num * edgetab_tmp = edgetab + verttab[vnodnbr];

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    const Vector<UInt> & connectivity = mesh.getConnectivity(type);

    for (UInt el = 0; el < nb_element; ++el) {
      *verttab_tmp = *(verttab_tmp - 1) + nb_nodes_per_element;
      verttab_tmp++;
      UInt * conn  = connectivity.values + el * nb_nodes_per_element;
      for (UInt i = 0; i < nb_nodes_per_element; ++i) {
	*(edgetab_tmp++) = *(conn++) + vnodbas;
      }
    }
  }

  SCOTCH_Mesh * meshptr = new SCOTCH_Mesh;

  SCOTCH_meshInit(meshptr);

  SCOTCH_meshBuild(meshptr,
		   velmbas, vnodbas,
		   velmnbr, vnodnbr,
		   verttab, vendtab,
		   velotab, vnlotab,
		   vlbltab,
		   edgenbr, edgetab);

  /// Check the mesh
  AKANTU_DEBUG_ASSERT(SCOTCH_meshCheck(meshptr) == 0,
		      "Scotch mesh is not consistent");

#ifndef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    /// save initial graph
    FILE *fmesh = fopen("ScotchMesh.msh", "w");
    SCOTCH_meshSave(meshptr, fmesh);
    fclose(fmesh);

    /// write geometry file
    std::ofstream fgeominit;
    fgeominit.open("GeomMeshFile.xyz");
    fgeominit << spatial_dimension << std::endl << nb_nodes << std::endl;

    const Vector<Real> & nodes = mesh.getNodes();
    Real * nodes_val = nodes.values;
    for (UInt i = 0; i < nb_nodes; ++i) {
      fgeominit << i << " ";
      for (UInt s = 0; s < spatial_dimension; ++s)
	fgeominit << *(nodes_val++) << " ";
      fgeominit << std::endl;;
    }
    fgeominit.close();
  }
#endif

  AKANTU_DEBUG_OUT();
  return meshptr;
}

/* -------------------------------------------------------------------------- */
void MeshPartitionScotch::destroyMesh(SCOTCH_Mesh * meshptr) {
  AKANTU_DEBUG_IN();

  SCOTCH_Num velmbas,  vnodbas,
             vnodnbr,  velmnbr,
            *verttab, *vendtab,
            *velotab, *vnlotab,
            *vlbltab,
             edgenbr, *edgetab,
             degrptr;


  SCOTCH_meshData(meshptr,
		  &velmbas, &vnodbas,
		  &velmnbr, &vnodnbr,
		  &verttab, &vendtab,
		  &velotab, &vnlotab,
		  &vlbltab,
		  &edgenbr, &edgetab,
		  &degrptr);

  delete [] verttab;
  delete [] edgetab;

  SCOTCH_meshExit(meshptr);

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
  SCOTCH_Num * edlotab = NULL;                  //integer  load  associated with
						//every edge ( optional )
  SCOTCH_Num * edgetab = dadjncy.values;        // adjacency array of every vertex
  SCOTCH_Num * vlbltab = NULL;                  // vertex label array (optional)

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

/* -------------------------------------------------------------------------- */
void MeshPartitionScotch::reorder() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Reordering the mesh " << mesh.getID());
  SCOTCH_Mesh * scotch_mesh = createMesh();


  destroyMesh(scotch_mesh);

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
