/**
 * @file   mesh_io_diana.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 *
 * @date   Sat Mar 26 20:43:38 2011
 *
 * @brief  handles diana meshes
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
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "mesh_io_diana.hh"

/* -------------------------------------------------------------------------- */
#include <string.h>
/* -------------------------------------------------------------------------- */
#include <stdio.h>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/*   Methods Implentations                                                    */
/* -------------------------------------------------------------------------- */

MeshIODiana::MeshIODiana() {
  canReadSurface      = true;
  canReadExtendedData = true;
  _diana_to_akantu_element_types["TE12L"] = _tetrahedron_4;
  _diana_to_akantu_element_types["HX24L"] = _hexahedron_8;
  _diana_to_akantu_mat_prop["YOUNG"] = "E";
  _diana_to_akantu_mat_prop["DENSIT"] = "rho";
  _diana_to_akantu_mat_prop["POISON"] = "nu";
}

/* -------------------------------------------------------------------------- */
MeshIODiana::~MeshIODiana() {
  std::map<std::string, Array<UInt> *>::iterator ng_it;
  std::map<std::string, std::vector<Element> *>::iterator eg_it;

  for (ng_it = node_groups.begin(); ng_it != node_groups.end(); ++ng_it) {
    delete ng_it->second;
  }

  for (eg_it = element_groups.begin(); eg_it != element_groups.end(); ++eg_it) {
    delete eg_it->second;
  }

}

/* -------------------------------------------------------------------------- */
inline void my_getline(std::ifstream & infile, std::string & line) {
  std::getline(infile, line); //read the line
  size_t pos = line.find("\r"); /// remove the extra \r if needed
  line = line.substr(0, pos);
}


/* -------------------------------------------------------------------------- */
void MeshIODiana::read(const std::string & filename, Mesh & mesh) {
  AKANTU_DEBUG_IN();

  std::ifstream infile;
  infile.open(filename.c_str());

  std::string line;
  UInt first_node_number = std::numeric_limits<UInt>::max();
  std::vector<Element> global_to_local_index;

  if(!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }

  while(infile.good()) {
    my_getline(infile, line);

    /// read all nodes
    if(line == "'COORDINATES'") {
      line = readCoordinates(infile, mesh, first_node_number);
    }

    /// read all elements
    if (line == "'ELEMENTS'") {
      line = readElements(infile, mesh, global_to_local_index, first_node_number);
    }

    /// read the material properties and write a .dat file
    if (line == "'MATERIALS'") {
      line = readMaterial(infile, filename);
    }

    /// read the material properties and write a .dat file
    if (line == "'GROUPS'") {
      line = readGroups(infile, global_to_local_index, first_node_number);
    }

  }
  infile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshIODiana::write(__attribute__((unused)) const std::string & filename,
			__attribute__((unused)) const Mesh & mesh) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readCoordinates(std::ifstream & infile, Mesh & mesh, UInt & first_node_number) {
  AKANTU_DEBUG_IN();

  Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());

  std::string line;

  UInt index;
  Real coord[3];

  do {
    my_getline(infile, line);
    if("'ELEMENTS'" == line)
      break;
    //end = true;
    //else {
    /// for each node, read the coordinates

    std::stringstream sstr_node(line);

    sstr_node >> index >> coord[0] >> coord[1] >> coord[2];

    //if (!sstr_node.fail())
    //break;

    first_node_number = first_node_number < index ? first_node_number : index;

    nodes.push_back(coord);
    // }
  } while(true);//!end);

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
UInt MeshIODiana::readInterval(std::stringstream & line,
			       std::set<UInt> & interval) {
  UInt first;
  line >> first;
  if(line.fail()) { return 0; }
  interval.insert(first);

  UInt second;
  int dash;
  dash = line.get();
  if(dash == '-') {
    line >> second;
    interval.insert(second);
    return 2;
  }

  if(line.fail())
    line.clear(std::ios::eofbit);  // in case of get at end of the line
  else line.unget();
  return 1;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readGroups(std::ifstream & infile,
				      std::vector<Element> & global_to_local_index,
				      UInt first_node_number) {
  AKANTU_DEBUG_IN();

  std::string line;
  my_getline(infile, line);

  bool reading_nodes_group = false;

  while(line != "'SUPPORTS'") {
    if(line == "NODES") {
      reading_nodes_group   = true;
      //      std::cout << line << std::endl;
      my_getline(infile, line);
    }

    if(line == "ELEMEN") {
      reading_nodes_group   = false;
      //      std::cout << line << std::endl;
      my_getline(infile, line);
    }

    std::stringstream *str = new std::stringstream(line);

    UInt id;
    std::string name;
    char c;
    *str >> id >> name >> c;
    //    std::cout << id << " " << name << " " << c << std::endl;

    Array<UInt> * list_ids = new Array<UInt>(0,1);

    UInt s = 1; bool end = false;
    while(!end) {
      while(!str->eof() && s != 0) {
	std::set<UInt> interval;
	s = readInterval(*str, interval);
	std::set<UInt>::iterator it = interval.begin();
	if(s == 1) list_ids->push_back(*it);
	if(s == 2) {
	  UInt first = *it;
	  ++it;
	  UInt second = *it;
	  for(UInt i = first; i <= second; ++i) {
	    list_ids->push_back(i);
	  }
	}
      }
      if(str->fail()) end = true;
      else {
	my_getline(infile, line);
	delete str;
	str = new std::stringstream(line);
      }
    }

    delete str;

    if(reading_nodes_group) {
      for (UInt i = 0; i < list_ids->getSize(); ++i) {
	(*list_ids)(i) -= first_node_number;
      }
      node_groups[name] = list_ids;
    } else {
      std::vector<Element> * elem = new std::vector<Element>;
      elem->reserve(list_ids->getSize());
      for (UInt i = 0; i < list_ids->getSize(); ++i) {
	elem->push_back(global_to_local_index[(*list_ids)(i)-1]);
      }

      element_groups[name] = elem;
      delete list_ids;
    }

    my_getline(infile, line);
  }

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readElements(std::ifstream & infile,
				      Mesh & mesh,
				      std::vector<Element> & global_to_local_index,
				      UInt first_node_number) {
  AKANTU_DEBUG_IN();

  std::string line;
  my_getline(infile, line);

  if("CONNECTIVITY" == line) {
    line = readConnectivity(infile, mesh, global_to_local_index, first_node_number);
  }

  /// read the line corresponding to the materials
  if ("MATERIALS" == line) {
    line = readMaterialElement(infile, mesh, global_to_local_index);
  }

  AKANTU_DEBUG_OUT();
  return line;
}


/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readConnectivity(std::ifstream & infile,
					  Mesh & mesh,
					  std::vector<Element> & global_to_local_index,
					  UInt first_node_number) {
  AKANTU_DEBUG_IN();

  Int index;
  std::string lline;

  std::string diana_type;
  ElementType akantu_type, akantu_type_old = _not_defined;
  Array<UInt> *connectivity = NULL;
  UInt node_per_element = 0;
  Element elem;
  UInt nb_elements_type = 0;

  bool end = false;
  do {
    my_getline(infile, lline);
    std::stringstream sstr_elem(lline);
    if(lline == "MATERIALS") end = true;
    else {
      /// traiter les coordonnees
      sstr_elem >> index;
      sstr_elem >> diana_type;

      akantu_type = _diana_to_akantu_element_types[diana_type];

      if(akantu_type != akantu_type_old) {
	connectivity = mesh.getConnectivityPointer(akantu_type);
	connectivity->resize(0);

	node_per_element = connectivity->getNbComponent();

	akantu_type_old = akantu_type;
	nb_elements_type = 0;
      }

      UInt local_connect[node_per_element];
      for(UInt j = 0; j < node_per_element; ++j) {
	UInt node_index;
	sstr_elem >> node_index;

	node_index -= first_node_number;
	local_connect[j] = node_index;
      }
      connectivity->push_back(local_connect);

      elem.type = akantu_type;
      elem.element = nb_elements_type;
      global_to_local_index.push_back(elem);
      nb_elements_type++;
    }
  }
  while(!end);

  AKANTU_DEBUG_OUT();
  return lline;
}

/* -------------------------------------------------------------------------- */
// UInt MeshIODiana::readInterval(std::stringstream & line,
// 			       std::set<UInt> & interval) {
//   UInt first;
//   line >> first;
//   if(line.fail()) { return 0; }
//   interval.insert(first);

//   UInt second;
//   char ignored;
//   line >> ignored >> second;
//   if(line.fail()) { line.clear(); return 1; }
//   interval.insert(second);
//   return 2;
// }

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readMaterialElement(std::ifstream & infile,
					     Mesh & mesh,
					     std::vector<Element> & global_to_local_index) {
  AKANTU_DEBUG_IN();


  // Array<UInt> vector_elements(nb_elements,1);
  //  ElementType akantu_type;
  std::string line;
  //  bool end = false;
  //  bool end_range = false;
  std::stringstream sstr_tag_name; sstr_tag_name << "tag_" << 0;
  //Array<UInt> * data = mesh.getUIntDataPointer(akantu_type, sstr_tag_name.str());

  Mesh::type_iterator it  = mesh.firstType();
  Mesh::type_iterator end = mesh.lastType();
  for(; it != end; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    mesh.getUIntDataPointer(*it, "material", _not_ghost)->resize(nb_element);
  }

  my_getline(infile, line);
  while(line != "'MATERIALS'") {
    line = line.substr(line.find('/') + 1, std::string::npos); // erase the first slash / of the line
    char tutu[250];
    strcpy(tutu, line.c_str());
    Array<UInt> temp_id(0, 2);
    UInt mat;
    while(true){
      std::stringstream sstr_intervals_elements(line);
      UInt id[2];
      char temp;
      while(sstr_intervals_elements.good()) {
	sstr_intervals_elements >> id[0] >> temp >> id[1]; // >> "/" >> mat;
	if(!sstr_intervals_elements.fail())
	  temp_id.push_back(id);
      }
      if (sstr_intervals_elements.fail()) {
	sstr_intervals_elements.clear();
	sstr_intervals_elements.ignore();
	sstr_intervals_elements >> mat;
	break;
      }
      my_getline(infile, line);
    }

    // loop over elements
    //    UInt * temp_id_val = temp_id.values;
    for (UInt i = 0; i < temp_id.getSize(); ++i)
      for (UInt j=temp_id(i,0); j<=temp_id(i,1); ++j) {
	Element & element = global_to_local_index[j - 1];
	UInt elem = element.element;
	ElementType type = element.type;
	Array<UInt> & data = *(mesh.getUIntDataPointer(type, "material", _not_ghost));
	data(elem) = mat;
      }

    my_getline(infile, line);
  }

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readMaterial(std::ifstream & infile,
				      const std::string & filename) {
  AKANTU_DEBUG_IN();

  std::stringstream mat_file_name;

  mat_file_name << "material_" << filename;

  std::ofstream material_file;
  material_file.open(mat_file_name.str().c_str());//mat_file_name.str());

  UInt mat_index;
  std::string line;

  bool first_mat = true;
  bool end = false;

  do{
    my_getline(infile, line);
    std::stringstream sstr_material(line);
    if("'GROUPS'" == line) {
      material_file << "]" << std::endl;
      end = true;
    }
    else {
      /// traiter les caractéristiques des matériaux
      sstr_material >> mat_index;

      if(!sstr_material.fail()) {
	if(!first_mat) {
	  material_file << "]" << std::endl;
	}
	material_file << "material elastic [" << std::endl;
	first_mat = false;
      } else {
	sstr_material.clear();
      }

      std::string prop_name;
      sstr_material >> prop_name;

      std::map<std::string, std::string>::iterator it;
      it = _diana_to_akantu_mat_prop.find(prop_name);

      if(it != _diana_to_akantu_mat_prop.end()) {
	Real value;
	sstr_material >> value;
	material_file << "\t" << it->second << " = " << value << std::endl;
      } else {
	//AKANTU_DEBUG_WARNING("In material reader, property " << it->first << "not recognized");
      }
    }
  } while (!end);

  AKANTU_DEBUG_OUT();
  return line;
}

__END_AKANTU__
