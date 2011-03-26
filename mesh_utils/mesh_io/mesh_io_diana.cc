/**
 * @file   mesh_io_diana.cc
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Mar 10 15:42:24 2011
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
#include <fstream>

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
    std::getline(infile, line);
    
    /// read all nodes
    if(line.find("COORDINATES") != std::string::npos) {
      readCoordinates(infile, mesh, first_node_number);
    }      
    
    /// read all elements
    if(line.find("CONNECTIVITY") != std::string::npos) {
      readConnectivity(infile, mesh, global_to_local_index, first_node_number);
    }
    
    /// read the line corresponding to the materials
    if (line.find("MATERIALS") != std::string::npos) {
      readMaterialElement(infile, mesh, global_to_local_index);
    }

    /// read the material properties and write a .dat file
    if (line.find("'MATERIALS'") != std::string::npos) {
      readMaterial(infile, mesh, filename);
    }
  }
  infile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshIODiana::write(const std::string & filename, const Mesh & mesh) {


}

/* -------------------------------------------------------------------------- */
void MeshIODiana::readCoordinates(std::ifstream & infile, Mesh & mesh, UInt & first_node_number) {
  AKANTU_DEBUG_IN();
  
  Vector<Real> & nodes = const_cast<Vector<Real> &>(mesh.getNodes());

  std::string line;

  UInt index;
  Real coord[3];
  bool end = false;

  do {
    std::getline(infile, line);
    //    if(line.find("'ELEMENTS'") != std::string::npos) end = true;
    //else {
      /// for each node, read the coordinates
      std::stringstream sstr_node(line);
      sstr_node >> index >> coord[0] >> coord[1] >> coord[2];	  
      
      first_node_number = first_node_number < index ? first_node_number : index;
  
      nodes.push_back(coord);
      // }
  } while(!sstr_node.fail());//!end);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshIODiana::readConnectivity(std::ifstream & infile, Mesh & mesh, std::vector<Element> & global_to_local_index, UInt first_node_number) {
  AKANTU_DEBUG_IN();

  Int index;
  std::string diana_type;
  ElementType akantu_type, akantu_type_old = _not_defined;
  Vector<UInt> *connectivity = NULL;
  UInt node_per_element = 0;
  Element elem;
  UInt nb_elements_type = 0;
   
  bool end = false;
  do {
    std::getline(infile, line);
    std::stringstream sstr_elem(line);
    if(line.find("MATERIALS") != std::string::npos) end = true;
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
}

/* -------------------------------------------------------------------------- */
void MeshIODiana::readMaterialElement(std::ifstream & infile, Mesh & mesh, std::vector<Element> & global_to_local_index) {
  AKANTU_DEBUG_IN();
  
  Vector<UInt> vector_elements(nb_elements,1);
  bool end = false;
  bool end_range = false;
  std::stringstream sstr_tag_name; sstr_tag_name << "tag_" << 0;
  Vector<UInt> * data = mesh.getUIntDataPointer(akantu_type, sstr_tag_name.str());
  
  do {
    std::getline(infile, line);
    if(line.find(" 'MATERIALS' ") != std::string::npos) 
      end = true;
    else {
      line.erase(0,0);
      /// erase the first slash / of the line
      Vector<UInt> temp_id;
      UInt mat;
      while(true){
	std::stringstream sstr_intervals_elements(line);
	UInt id[2];
	while(sstr_intervals_elements.good()) {
	  sstr_intervals_elements >> id[0] >> "-" >> id[1]; // >> "/" >> mat;
	  temp_id.push_back(id);
	}
	if (sstr_intervals_elements.fail()) {
	  sstr_intervals_elements >> mat;
	  break;
	}
	std::getline(infile, line);
      }
      /// loop over elements
      UInt * temp_id_val = temp_id.values;
      for (UInt i = 0; i < temp_id.size()/2; ++i)
	for (UInt j=temp_id_val[i]; j=temp_id_val[i+1];j+=2)
	  mesh.getUIntDataPointer(global_to_local_index[j - 1].type, "material")->at(global_to_local_index[j - 1].element) = mat;
    }
  }
  while(!end);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshIODiana::readMaterial(std::ifstream & infile, Mesh & mesh, const std::string & filename) {
  AKANTU_DEBUG_IN();
  
  std::stringstream mat_file_name;
  
  mat_file_name << "material_" << filename;
  
  std::ofstream material_file;
  material_file.open(mat_file_name.str());

  UInt mat_index;
  std::string line;
  
  bool first_mat = true;
  bool end = false;
  
  do{
    std::getline(infile, line);
    std::stringstream sstr_material(line);
    if(line.find("GROUPS") != std::string::npos) {
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
	AKANTU_DEBUG_WARNING("In material reader, property " << it->first << "not recognized");
      }
    }
  } while (!end);
  
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

