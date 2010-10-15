/**
 * @file   mesh_io_msh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 18 11:36:35 2010
 *
 * @brief  Read/Write for MSH files generated by gmsh
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -----------------------------------------------------------------------------
   Version (Legacy) 1.0

   $NOD
   number-of-nodes
   node-number x-coord y-coord z-coord
   ...
   $ENDNOD
   $ELM
   number-of-elements
   elm-number elm-type reg-phys reg-elem number-of-nodes node-number-list
   ...
   $ENDELM
   -----------------------------------------------------------------------------
   Version 2.1

   $MeshFormat
   version-number file-type data-size
   $EndMeshFormat
   $Nodes
   number-of-nodes
   node-number x-coord y-coord z-coord
   ...
   $EndNodes
   $Elements
   number-of-elements
   elm-number elm-type number-of-tags < tag > ... node-number-list
   ...
   $EndElements
   $PhysicalNames
   number-of-names
   physical-dimension physical-number "physical-name"
   ...
   $EndPhysicalNames
   $NodeData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   node-number value ...
   ...
   $EndNodeData
   $ElementData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   elm-number value ...
   ...
   $EndElementData
   $ElementNodeData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   elm-number number-of-nodes-per-element value ...
   ...
   $ElementEndNodeData

   -----------------------------------------------------------------------------
   elem-type

   1:  2-node line.
   2:  3-node triangle.
   3:  4-node quadrangle.
   4:  4-node tetrahedron.
   5:  8-node hexahedron.
   6:  6-node prism.
   7:  5-node pyramid.
   8:  3-node second order line
   9:  6-node second order triangle
   10: 9-node second order quadrangle
   11: 10-node second order tetrahedron
   12: 27-node second order hexahedron
   13: 18-node second order prism
   14: 14-node second order pyramid
   15: 1-node point.
   16: 8-node second order quadrangle
   17: 20-node second order hexahedron
   18: 15-node second order prism
   19: 13-node second order pyramid
   20: 9-node third order incomplete triangle
   21: 10-node third order triangle
   22: 12-node fourth order incomplete triangle
   23: 15-node fourth order triangle
   24: 15-node fifth order incomplete triangle
   25: 21-node fifth order complete triangle
   26: 4-node third order edge
   27: 5-node fourth order edge
   28: 6-node fifth order edge
   29: 20-node third order tetrahedron
   30: 35-node fourth order tetrahedron
   31: 56-node fifth order tetrahedron
   -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "mesh_io_msh.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/*   Constant arrays initialisation                                           */
/* -------------------------------------------------------------------------- */
UInt MeshIOMSH::_read_order[_max_element_type][MAX_NUMBER_OF_NODE_PER_ELEMENT] = {
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // _not_defined
  { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, // _line1
  { 0, 1, 2, 0, 0, 0, 0, 0, 0, 0 }, // _line2
  { 0, 1, 2, 0, 0, 0, 0, 0, 0, 0 }, // _triangle_1
  { 0, 1, 2, 3, 4, 5, 6, 0, 0, 0 }, // _triangle_2
  { 0, 1, 2, 3, 4, 0, 0, 0, 0, 0 }, // _tetrahedra_1
  { 0, 1, 2, 3, 4, 5, 6, 7, 9, 8 }, // _tetrahedra_2
};

UInt MeshIOMSH::_msh_nodes_per_elem[16] =
  { 0, // element types began at 1
    2, 3, 4, 4,  8,  6,  5,  // 1st order
    3, 6, 9, 10, 27, 18, 14, // 2st order
    1
  };

ElementType MeshIOMSH::_msh_to_akantu_element_types[16] =
  { _not_defined, // element types began at 1
    _line_1,      _triangle_1,  _not_defined, _tetrahedra_1, // 1st order
    _not_defined, _not_defined, _not_defined,
    _line_2,      _triangle_2,  _not_defined, _tetrahedra_2, // 2nd order
    _not_defined, _not_defined, _not_defined,
    _not_defined
  };

MeshIOMSH::MSHElementType MeshIOMSH::_akantu_to_msh_element_types[_max_element_type
] =
  { _msh_not_defined,   // _not_defined
    _msh_line_1,        // _line_1
    _msh_line_2,        // _line_2
    _msh_triangle_1,    // _triangle_1
    _msh_triangle_2,    // _triangle_2
    _msh_tetrahedron_1, // _tetrahedra_1
    _msh_tetrahedron_2  // _tetrahedra_2
  };


/* -------------------------------------------------------------------------- */
/*   Methods Implentations                                                    */
/* -------------------------------------------------------------------------- */

MeshIOMSH::MeshIOMSH() {
  canReadSurface      = false;
  canReadExtendedData = false;
}

/* -------------------------------------------------------------------------- */
MeshIOMSH::~MeshIOMSH() {

}

/* -------------------------------------------------------------------------- */
void MeshIOMSH::read(const std::string & filename, Mesh & mesh) {

  std::ifstream infile;
  infile.open(filename.c_str());

  std::string line;
  UInt first_node_number = 0, last_node_number = 0,
    file_format = 1, current_line = 0;


  if(!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }

  while(infile.good()) {
    std::getline(infile, line);
    current_line++;

    /// read the header
    if(line == "$MeshFormat") {
      std::getline(infile, line); /// the format line
      std::getline(infile, line); /// the end of block line
      current_line += 2;
      file_format = 2;
    }

    /// read all nodes
    if(line == "$Nodes" || line == "$NOD") {
      UInt nb_nodes;

      std::getline(infile, line);
      std::stringstream sstr(line);
      sstr >> nb_nodes;
      current_line++;

      Vector<Real> & nodes = const_cast<Vector<Real> &>(mesh.getNodes());
      nodes.resize(nb_nodes);

      UInt index;
      Real coord[3];
      UInt spatial_dimension = nodes.getNbComponent();
      /// for each node, read the coordinates
      for(UInt i = 0; i < nb_nodes; ++i) {
	UInt offset = i * spatial_dimension;

	std::getline(infile, line);
	std::stringstream sstr_node(line);
	sstr_node >> index >> coord[0] >> coord[1] >> coord[2];
	current_line++;

	first_node_number = first_node_number < index ? first_node_number : index;
	last_node_number  = last_node_number  > index ? last_node_number  : index;

	/// read the coordinates
	for(UInt j = 0; j < spatial_dimension; ++j)
	  nodes.values[offset + j] = coord[j];
      }
      std::getline(infile, line); /// the end of block line
    }


    /// read all elements
    if(line == "$Elements" || line == "$ELM") {
      UInt nb_elements;

      std::getline(infile, line);
      std::stringstream sstr(line);
      sstr >> nb_elements;
      current_line++;

      Int index;
      UInt msh_type;
      ElementType akantu_type, akantu_type_old = _not_defined;
      Vector<UInt> *connectivity = NULL;
      UInt nb_elements_read = 0;
      UInt node_per_element = 0;

      for(UInt i = 0; i < nb_elements; ++i) {
	std::getline(infile, line);
	std::stringstream sstr_elem(line);
	current_line++;

	sstr_elem >> index;
	sstr_elem >> msh_type;

	/// get the connectivity vector depending on the element type
	akantu_type = _msh_to_akantu_element_types[msh_type];

	if(akantu_type == _not_defined) continue;

	if(akantu_type != akantu_type_old) {
	  if(connectivity)
	    connectivity->resize(nb_elements_read);

	  connectivity = mesh.getConnectivityPointer(akantu_type);
	  connectivity->resize(nb_elements);

	  node_per_element = connectivity->getNbComponent();
	  akantu_type_old = akantu_type;
	  nb_elements_read = 0;
	}

	/// read tags informations
	if(file_format == 2) {
	  UInt nb_tags;
	  sstr_elem >> nb_tags;
	  for(UInt j = 0; j < nb_tags; ++j) {
	    Int tag;
	    sstr_elem >> tag; ///@todo read to get extended information on elements
	  }
	} else if (file_format == 1) {
	  Int tag;
	  sstr_elem >> tag; //reg-phys
	  sstr_elem >> tag; //reg-elem
	  sstr_elem >> tag; //number-of-nodes
	}

	/// read the connectivities informations
	UInt offset = nb_elements_read * node_per_element;
	for(UInt j = 0; j < node_per_element; ++j) {
	  UInt node_index;
	  sstr_elem >> node_index;

	  AKANTU_DEBUG_ASSERT(node_index <= last_node_number,
			     "Node number not in range : line " << current_line);

	  node_index -= first_node_number + 1;
	  connectivity->values[offset + _read_order[akantu_type][j]] = node_index;
	}
	nb_elements_read++;
      }
      connectivity->resize(nb_elements_read);
      std::getline(infile, line); /// the end of block line
    }

    if((line[0] == '$') && (line.find("End") == std::string::npos)) {
      AKANTU_DEBUG_WARNING("Unsuported block_kind " << line
			  << " at line " << current_line);
    }
  }

  infile.close();
}

/* -------------------------------------------------------------------------- */
void MeshIOMSH::write(const std::string & filename, const Mesh & mesh) {
  std::ofstream outfile;
  const Vector<Real> & nodes = mesh.getNodes();

  outfile.open(filename.c_str());

  outfile << "$MeshFormat" << std::endl;
  outfile << "2 0 8" << std::endl;;
  outfile << "$EndMeshFormat" << std::endl;;

  outfile << "$Nodes" << std::endl;;
  outfile << nodes.getSize() << std::endl;;

  outfile << std::uppercase;
  for(UInt i = 0; i < nodes.getSize(); ++i) {
    Int offset = i * nodes.getNbComponent();
    outfile << i+1;
    for(UInt j = 0; j < nodes.getNbComponent(); ++j) {
      outfile << " " << nodes.values[offset + j];
    }

    if(nodes.getNbComponent() == 2)
      outfile << " " << 0.;
    outfile << std::endl;;
  }
  outfile << std::nouppercase;
  outfile << "$EndNodes" << std::endl;;


  outfile << "$Elements" << std::endl;;

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  Int nb_elements = 0;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    const Vector<UInt> & connectivity = mesh.getConnectivity(*it);
    nb_elements += connectivity.getSize();
  }
  outfile << nb_elements << std::endl;

  UInt element_idx = 1;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    const Vector<UInt> & connectivity = mesh.getConnectivity(type);

    for(UInt i = 0; i < connectivity.getSize(); ++i) {
      UInt offset = i * connectivity.getNbComponent();
      outfile << element_idx << " " << type << " 3 0 0 0";

      for(UInt j = 0; j < connectivity.getNbComponent(); ++j) {
	outfile << " " << connectivity.values[offset + j];
      }
      outfile << std::endl;
    }
    element_idx++;
  }

  outfile << "$EndElements" << std::endl;;

  outfile.close();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__

