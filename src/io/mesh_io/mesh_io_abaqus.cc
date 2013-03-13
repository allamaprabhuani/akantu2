/**
 * @file   mesh_io_abaqus.cc
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Fri Oct  7 10:58:00 2011
 *
 * @brief  read a mesh from an abaqus input file
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

// std library header files
#include <fstream>

// akantu header files
#include "mesh_io_abaqus.hh"
#include "aka_ci_string.hh"
#include "aka_abaqus_parser.hh"

__BEGIN_AKANTU__

void MeshIOAbaqus::read(const std::string& filename, Mesh& mesh) {
    
    
    typedef ci_string string_type;
    typedef const string_type& key_type;
    typedef string_type::value_type char_type;
    typedef std::basic_ifstream<char_type, ci_char_traits> stream_type;
    typedef std::basic_istringstream<char_type, ci_char_traits> sstream_type;
    typedef std::list<string_type> token_container;
    typedef token_container::iterator token_iterator;
    
    // problem dimension
    UInt dim = 0;
    
    // input stream object
    stream_type ifs;
    
    ifs.open(filename.c_str(), std::ifstream::in);
    
    if (!ifs)
        throw debug::Exception("Could not open file", __FILE__, __LINE__);
    
    string_type line;
    token_container tokens;
    unsigned long lineno = 0;
    
    // start loop over lines in the file
    while (ifs.good()) {
        
        // process line, eliminating comment lines and empty lines
        if (!getline(ifs, line, lineno))
            continue;
        
        // get keyword
        boost::char_separator<char_type> sep("*,\r\n");
        string_type keyword = first_token(line, sep);
        
        // node section
        if (keyword == "node") {
            
            Array<Real>& nodes = const_cast<Array<Real> &>(mesh.getNodes());
            
            // read nodes from file
            while (ifs.peek() != '*') {
                
                // process line
                if (!getline(ifs, line, lineno))
                    continue;
                
                // tokenize line
                tokenize(line, tokens);
                
                Real coord[] = {0., 0., 0.};
                
                int i = 0;
                for (token_iterator it = ++tokens.begin(); it != tokens.end();)
                    coord[i++] = stream_cast<Real>(*it++, lineno);
                
                // push coordinate
                nodes.push_back(coord);
            }      
        } // node section
        
        // node section
        else if (keyword == "element") {
            
            // tokenize line
            tokenize(line, tokens);
            
            // get value for key
            key_type eltype = value_for_key("TYPE", tokens);
            
            // get element dimensionality
            Abaqus_element_map abaqus_map;
            
            
            // read nodes from file
            while (ifs.peek() != '*') {
                
                // process line
                if (!getline(ifs, line, lineno))
                    continue;
                
                // STOPPED HERE FOR NOW
                assert(false);
                
                //        Int index;
                //        std::string lline;
                //        
                //        std::string diana_type;
                //        ElementType akantu_type, akantu_type_old = _not_defined;
                //        Array<UInt> *connectivity = NULL;
                //        UInt node_per_element = 0;
                //        Element elem;
                //        UInt nb_elements_type = 0;
                //        
                //        bool end = false;
                //        do {
                //          my_getline(infile, lline);
                //          std::stringstream sstr_elem(lline);
                //          if(lline == "MATERIALS") end = true;
                //          else {
                //            /// traiter les coordonnees
                //            sstr_elem >> index;
                //            sstr_elem >> diana_type;
                //            
                //            akantu_type = _diana_to_akantu_element_types[diana_type];
                //            
                //            if(akantu_type != akantu_type_old) {
                //              connectivity = mesh.getConnectivityPointer(akantu_type);
                //              connectivity->resize(0);
                //              
                //              node_per_element = connectivity->getNbComponent();
                //              
                //              akantu_type_old = akantu_type;
                //              nb_elements_type = 0;
                //            }
                //            
                //            UInt local_connect[node_per_element];
                //            for(UInt j = 0; j < node_per_element; ++j) {
                //              UInt node_index;
                //              sstr_elem >> node_index;
                //              
                //              node_index -= first_node_number;
                //              local_connect[j] = node_index;
                //            }
                //            connectivity->push_back(local_connect);
                //            
                //            elem.type = akantu_type;
                //            elem.element = nb_elements_type;
                //            global_to_local_index.push_back(elem);
                //            nb_elements_type++;
                //          }
                //        }
                //        while(!end);
                
            }
            
        } // element section
        
        else {
            
            cout<<"Not implemented keyword "<<keyword<<endl;
            
        }
    }
}


__END_AKANTU__
