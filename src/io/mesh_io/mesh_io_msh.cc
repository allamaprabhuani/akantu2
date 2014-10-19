/**
 * @file   mesh_io_msh.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jul 04 2014
 *
 * @brief  Read/Write for MSH files generated by gmsh
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_io.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/*   Methods Implentations                                                    */
/* -------------------------------------------------------------------------- */

MeshIOMSH::MeshIOMSH() {
  canReadSurface      = true;
  canReadExtendedData = true;

  _msh_nodes_per_elem[_msh_not_defined   ] = 0;
  _msh_nodes_per_elem[_msh_segment_2     ] = 2;
  _msh_nodes_per_elem[_msh_triangle_3    ] = 3;
  _msh_nodes_per_elem[_msh_quadrangle_4  ] = 4;
  _msh_nodes_per_elem[_msh_tetrahedron_4 ] = 4;
  _msh_nodes_per_elem[_msh_hexahedron_8  ] = 8;
  _msh_nodes_per_elem[_msh_prism_1       ] = 6;
  _msh_nodes_per_elem[_msh_pyramid_1     ] = 1;
  _msh_nodes_per_elem[_msh_segment_3     ] = 3;
  _msh_nodes_per_elem[_msh_triangle_6    ] = 6;
  _msh_nodes_per_elem[_msh_quadrangle_9  ] = 9;
  _msh_nodes_per_elem[_msh_tetrahedron_10] = 10;
  _msh_nodes_per_elem[_msh_hexahedron_27 ] = 27;
  _msh_nodes_per_elem[_msh_prism_18      ] = 18;
  _msh_nodes_per_elem[_msh_pyramid_14    ] = 14;
  _msh_nodes_per_elem[_msh_point         ] = 1;
  _msh_nodes_per_elem[_msh_quadrangle_8  ] = 8;

  _msh_to_akantu_element_types[_msh_not_defined   ] = _not_defined;
  _msh_to_akantu_element_types[_msh_segment_2     ] = _segment_2;
  _msh_to_akantu_element_types[_msh_triangle_3    ] = _triangle_3;
  _msh_to_akantu_element_types[_msh_quadrangle_4  ] = _quadrangle_4;
  _msh_to_akantu_element_types[_msh_tetrahedron_4 ] = _tetrahedron_4;
  _msh_to_akantu_element_types[_msh_hexahedron_8  ] = _hexahedron_8;
  _msh_to_akantu_element_types[_msh_prism_1       ] = _pentahedron_6;
  _msh_to_akantu_element_types[_msh_pyramid_1     ] = _not_defined;
  _msh_to_akantu_element_types[_msh_segment_3     ] = _segment_3;
  _msh_to_akantu_element_types[_msh_triangle_6    ] = _triangle_6;
  _msh_to_akantu_element_types[_msh_quadrangle_9  ] = _not_defined;
  _msh_to_akantu_element_types[_msh_tetrahedron_10] = _tetrahedron_10;
  _msh_to_akantu_element_types[_msh_hexahedron_27 ] = _not_defined;
  _msh_to_akantu_element_types[_msh_prism_18      ] = _not_defined;
  _msh_to_akantu_element_types[_msh_pyramid_14    ] = _not_defined;
  _msh_to_akantu_element_types[_msh_point         ] = _point_1;
  _msh_to_akantu_element_types[_msh_quadrangle_8  ] = _quadrangle_8;

  _akantu_to_msh_element_types[_not_defined     ] = _msh_not_defined;
  _akantu_to_msh_element_types[_segment_2       ] = _msh_segment_2;
  _akantu_to_msh_element_types[_segment_3       ] = _msh_segment_3;
  _akantu_to_msh_element_types[_triangle_3      ] = _msh_triangle_3;
  _akantu_to_msh_element_types[_triangle_6      ] = _msh_triangle_6;
  _akantu_to_msh_element_types[_tetrahedron_4   ] = _msh_tetrahedron_4;
  _akantu_to_msh_element_types[_tetrahedron_10  ] = _msh_tetrahedron_10;
  _akantu_to_msh_element_types[_quadrangle_4    ] = _msh_quadrangle_4;
  _akantu_to_msh_element_types[_quadrangle_8    ] = _msh_quadrangle_8;
  _akantu_to_msh_element_types[_hexahedron_8    ] = _msh_hexahedron_8;
  _akantu_to_msh_element_types[_pentahedron_6   ] = _msh_prism_1;
  _akantu_to_msh_element_types[_point_1         ] = _msh_point;
#if defined(AKANTU_STRUCTURAL_MECHANICS)
  _akantu_to_msh_element_types[_bernoulli_beam_2] = _msh_segment_2;
  _akantu_to_msh_element_types[_bernoulli_beam_3] = _msh_segment_2;
  _akantu_to_msh_element_types[_kirchhoff_shell] = _msh_triangle_3;
#endif

  std::map<ElementType, MSHElementType>::iterator it;
  for(it = _akantu_to_msh_element_types.begin();
      it != _akantu_to_msh_element_types.end(); ++it) {
    UInt nb_nodes = _msh_nodes_per_elem[it->second];

    UInt * tmp = new UInt[nb_nodes];
    for (UInt i = 0; i < nb_nodes; ++i) {
      tmp[i] = i;
    }

    switch(it->first) {
    case _tetrahedron_10:
      tmp[8] = 9;
      tmp[9] = 8;
      break;
    default:
      //nothing to change
      break;
    }
    _read_order[it->first] = tmp;
  }
}

/* -------------------------------------------------------------------------- */
MeshIOMSH::~MeshIOMSH() {
  std::map<ElementType, MSHElementType>::iterator it;
  for(it = _akantu_to_msh_element_types.begin();
      it != _akantu_to_msh_element_types.end(); ++it) {
    delete [] _read_order[it->first];
  }
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

    /// read the physical names
    if(line == "$PhysicalNames") {
      std::getline(infile, line); /// the format line
      std::stringstream sstr(line);

      UInt num_of_phys_names;
      sstr >> num_of_phys_names;


      for(UInt k(0); k < num_of_phys_names; k++) {
        std::getline(infile, line);
        std::stringstream sstr_phys_name(line);
        UInt phys_name_id;
        UInt phys_dim;

        sstr_phys_name >> phys_dim >> phys_name_id;

	std::size_t b = line.find("\"");
	std::size_t e = line.rfind("\"");
        std::string phys_name = line.substr(b + 1, e - b -1);

        phys_name_map[phys_name_id] = phys_name;
      }
    }

    /// read all nodes
    if(line == "$Nodes" || line == "$NOD") {
      UInt nb_nodes;

      std::getline(infile, line);
      std::stringstream sstr(line);
      sstr >> nb_nodes;
      current_line++;

      Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());
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
          nodes.storage()[offset + j] = coord[j];
      }
      std::getline(infile, line); /// the end of block line
    }


    /// read all elements
    if(line == "$Elements" || line == "$ELM") {
      UInt nb_elements;

      UInt * read_order = NULL;

      std::getline(infile, line);
      std::stringstream sstr(line);
      sstr >> nb_elements;
      current_line++;

      Int index;
      UInt msh_type;
      ElementType akantu_type, akantu_type_old = _not_defined;
      Array<UInt> *connectivity = NULL;
      UInt node_per_element = 0;

      for(UInt i = 0; i < nb_elements; ++i) {
        std::getline(infile, line);
        std::stringstream sstr_elem(line);
        current_line++;

        sstr_elem >> index;
        sstr_elem >> msh_type;

        /// get the connectivity vector depending on the element type
        akantu_type = _msh_to_akantu_element_types[(MSHElementType) msh_type];

        if(akantu_type == _not_defined) {
          AKANTU_DEBUG_WARNING("Unsuported element kind " << msh_type
                               << " at line " << current_line);
          continue;
        }

        if(akantu_type != akantu_type_old) {
          connectivity = mesh.getConnectivityPointer(akantu_type);
//          connectivity->resize(0);

          node_per_element = connectivity->getNbComponent();
          akantu_type_old = akantu_type;
          read_order = _read_order[akantu_type];
        }

        /// read tags informations
        if(file_format == 2) {
          UInt nb_tags;
          sstr_elem >> nb_tags;
          for(UInt j = 0; j < nb_tags; ++j) {
            Int tag;
            sstr_elem >> tag;
            std::stringstream sstr_tag_name; sstr_tag_name << "tag_" << j;
            Array<UInt> * data = mesh.getDataPointer<UInt>(sstr_tag_name.str(), akantu_type, _not_ghost);
            data->push_back(tag);
          }
        } else if (file_format == 1) {
          Int tag;
          sstr_elem >> tag; //reg-phys
          std::string tag_name = "tag_0";
          Array<UInt> * data = mesh.getDataPointer<UInt>(tag_name, akantu_type, _not_ghost);
          data->push_back(tag);

          sstr_elem >> tag; //reg-elem
          tag_name = "tag_1";
          data = mesh.getDataPointer<UInt>(tag_name, akantu_type, _not_ghost);
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
          local_connect[read_order[j]] = node_index;
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

  mesh.updateTypesOffsets(_not_ghost);

  infile.close();

  if(!phys_name_map.empty()) {
    for(Mesh::type_iterator type_it = mesh.firstType(); type_it != mesh.lastType(); ++type_it) {
      Array<std::string> * name_vec = mesh.getDataPointer<std::string>("physical_names", *type_it);
      const Array<UInt> & tags_vec = mesh.getData<UInt>("tag_0", *type_it);

      for(UInt i(0); i < tags_vec.getSize(); i++) {
        std::map<UInt, std::string>::const_iterator map_it = phys_name_map.find(tags_vec(i));

        if(map_it == phys_name_map.end()) {
          std::stringstream sstm;
          sstm << tags_vec(i);
          name_vec->operator()(i) = sstm.str();
        } else {
          name_vec->operator()(i) = map_it->second;
        }
      }
    }
  }

  MeshUtils::fillElementToSubElementsData(mesh);
}

/* -------------------------------------------------------------------------- */
void MeshIOMSH::write(const std::string & filename, const Mesh & mesh) {
  std::ofstream outfile;
  const Array<Real> & nodes = mesh.getNodes();

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
      outfile << " " << nodes.storage()[offset + j];
    }

    for (UInt p = nodes.getNbComponent(); p < 3; ++p)
      outfile << " " << 0.;
    outfile << std::endl;;
  }
  outfile << std::nouppercase;
  outfile << "$EndNodes" << std::endl;;


  outfile << "$Elements" << std::endl;;

  Mesh::type_iterator it  = mesh.firstType(_all_dimensions, _not_ghost, _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType(_all_dimensions, _not_ghost, _ek_not_defined);

  Int nb_elements = 0;
  for(; it != end; ++it) {
    const Array<UInt> & connectivity = mesh.getConnectivity(*it, _not_ghost);
    nb_elements += connectivity.getSize();
  }
  outfile << nb_elements << std::endl;

  UInt element_idx = 1;
  for(it  = mesh.firstType(_all_dimensions, _not_ghost, _ek_not_defined); it != end; ++it) {
    ElementType type = *it;
    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    UInt * tag[2] = {NULL, NULL};
    try {
      const Array<UInt> & data_tag_0 = mesh.getData<UInt>("tag_0", type, _not_ghost);
      tag[0] = data_tag_0.storage();
    } catch(...) { tag[0] = NULL; }

    try {
      const Array<UInt> & data_tag_1 = mesh.getData<UInt>("tag_1", type, _not_ghost);
      tag[1] = data_tag_1.storage();
    } catch(...) { tag[1] = NULL; }

    for(UInt i = 0; i < connectivity.getSize(); ++i) {
      UInt offset = i * connectivity.getNbComponent();
      outfile << element_idx << " " << _akantu_to_msh_element_types[type] << " 2";

      /// \todo write the real data in the file
      for (UInt t = 0; t < 2; ++t)
        if(tag[t]) outfile << " " << tag[t][i];
        else outfile << " 0";

      for(UInt j = 0; j < connectivity.getNbComponent(); ++j) {
        outfile << " " << connectivity.storage()[offset + j] + 1;
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
