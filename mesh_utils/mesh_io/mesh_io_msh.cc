/**
 * @file   mesh_io_msh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 18 11:36:35 2010
 *
 * @brief  Read/Write for MSH files generated by gmsh
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
  { 0, 1, 2, 0, 0, 0, 0, 0, 0, 0 }, // _triangle_3
  { 0, 1, 2, 3, 4, 5, 6, 0, 0, 0 }, // _triangle_6
  { 0, 1, 2, 3, 4, 0, 0, 0, 0, 0 }, // _tetrahedron_4
  { 0, 1, 2, 3, 4, 5, 6, 7, 9, 8 }, // _tetrahedron_10
  { 0, 1, 2, 3, 0, 0, 0, 0, 0, 0 }, // _quadrangle_4
  { 0, 1, 2, 3, 4, 5, 6, 7, 0, 0 }, // _quadrangle_8
  { 0, 1, 2, 3, 4, 5, 6, 7, 0, 0 }  // _hexahedron_8
};

UInt MeshIOMSH::_msh_nodes_per_elem[17] =
  { 0, // element types began at 1
    2, 3, 4, 4,  8,  6,  5,  // 1st order
    3, 6, 9, 10, 27, 18, 14, // 2st order
    1, 8
  };

ElementType MeshIOMSH::_msh_to_akantu_element_types[17] =
  { _not_defined, // element types began at 1
    _segment_2,   _triangle_3,  _quadrangle_4, _tetrahedron_4,  // 1st order
    _hexahedron_8, _not_defined, _not_defined,
    _segment_3,   _triangle_6,  _not_defined,  _tetrahedron_10, // 2nd order
    _not_defined, _not_defined, _not_defined,
    _not_defined, _quadrangle_8
  };

MeshIOMSH::MSHElementType MeshIOMSH::_akantu_to_msh_element_types[_max_element_type] =
  {
    _msh_not_defined,   // _not_defined
    _msh_segment_2,     // _segment_2
    _msh_segment_3,     // _segment_3
    _msh_triangle_3,    // _triangle_3
    _msh_triangle_6,    // _triangle_6
    _msh_tetrahedron_4, // _tetrahedron_4
    _msh_tetrahedron_10,// _tetrahedron_10
    _msh_quadrangle_4,  // _quadrangle_4
    _msh_quadrangle_8,  // _quadrangle_4
    _msh_hexahedron_8   // _hexahedron_8
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
  UInt first_node_number = std::numeric_limits<UInt>::max(), last_node_number = 0,
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
      std::stringstream sstr(line);
      std::string version; sstr >> version;
      Int format; sstr >> format;
      if(format != 0) AKANTU_DEBUG_ERROR("This reader can only read ASCII files.");
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
      mesh.nb_global_nodes = nb_nodes;


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

	first_node_number = std::min(first_node_number,index);
	last_node_number  = std::max(last_node_number, index);

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
	  connectivity = mesh.getConnectivityPointer(akantu_type);
	  connectivity->resize(0);

	  node_per_element = connectivity->getNbComponent();
	  akantu_type_old = akantu_type;
	}

	/// read tags informations
	if(file_format == 2) {
	  UInt nb_tags;
	  sstr_elem >> nb_tags;
	  for(UInt j = 0; j < nb_tags; ++j) {
	    Int tag;
	    sstr_elem >> tag;
	    std::stringstream sstr_tag_name; sstr_tag_name << "tag_" << j;
	    Vector<UInt> * data = mesh.getUIntDataPointer(akantu_type, sstr_tag_name.str(), _not_ghost);
	    data->push_back(tag);
	  }
	} else if (file_format == 1) {
	  Int tag;
	  sstr_elem >> tag; //reg-phys
	  std::string tag_name = "tag_0";
	  Vector<UInt> * data = mesh.getUIntDataPointer(akantu_type, tag_name, _not_ghost);
	  data->push_back(tag);

	  sstr_elem >> tag; //reg-elem
	  tag_name = "tag_1";
	  data = mesh.getUIntDataPointer(akantu_type, tag_name, _not_ghost);
	  data->push_back(tag);

	  sstr_elem >> tag; //number-of-nodes
	}

	UInt local_connect[node_per_element];
	for(UInt j = 0; j < node_per_element; ++j) {
	  UInt node_index;
	  sstr_elem >> node_index;

	  AKANTU_DEBUG_ASSERT(node_index <= last_node_number,
			     "Node number not in range : line " << current_line);

	  node_index -= first_node_number;
	  local_connect[_read_order[akantu_type][j]] = node_index;
	}
	connectivity->push_back(local_connect);
      }
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
  outfile << "2.1 0 8" << std::endl;;
  outfile << "$EndMeshFormat" << std::endl;;


  outfile << std::setprecision(std::numeric_limits<Real>::digits10);
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
    const Vector<UInt> & connectivity = mesh.getConnectivity(*it, _not_ghost);
    nb_elements += connectivity.getSize();
  }
  outfile << nb_elements << std::endl;

  UInt element_idx = 1;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    const Vector<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);


    for(UInt i = 0; i < connectivity.getSize(); ++i) {
      UInt offset = i * connectivity.getNbComponent();
      outfile << element_idx << " " << _akantu_to_msh_element_types[type] << " 3 1 1 0"; /// \todo write the real data in the file

      for(UInt j = 0; j < connectivity.getNbComponent(); ++j) {
	outfile << " " << connectivity.values[offset + j] + 1;
      }
      outfile << std::endl;
      element_idx++;
    }
  }

  outfile << "$EndElements" << std::endl;;

  outfile.close();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__

