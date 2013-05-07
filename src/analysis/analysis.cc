/**
 * @file   analysis.cc
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Wed Oct  5 15:35:00 2011
 *
 * @brief  analysis implementation
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
#include <iostream>
#include <fstream>


// akantu header files
#include "analysis.hh"
#include "aka_error.hh"
#include "aka_ci_string.hh"
#include "aka_abaqus_parser.hh"


__BEGIN_AKANTU__

using std::cout;
using std::endl;


void Analysis<Abaqus>::read_file(const char* filename) {
  
  typedef ci_string string_type;
  typedef const string_type& key_type;
  typedef string_type::value_type char_type;
  typedef std::basic_ifstream<char_type, ci_char_traits> stream_type;
  typedef std::basic_istringstream<char_type, ci_char_traits> sstream_type;
  typedef std::list<string_type> token_container;
  typedef token_container::iterator token_iterator;
  
  
  // problem dimension
//  UInt dim = 0;
  
  // input stream object
  stream_type ifs;
  
  ifs.open(filename, std::ifstream::in );
  
  if (!ifs)
    throw debug::Exception("Could not open file", __FILE__, __LINE__);
  
  string_type line;
  token_container tokens;
  unsigned long lineno = 0;
    
  bool mesh_created = false;
  
  // start loop over lines in the file
  while (ifs.good()) {
    
    // process line, eliminating comment lines and empty lines
    if (!getline(ifs, line, lineno))
      continue;
    
    // get keyword
    boost::char_separator<char_type> sep("*,\r\n");
    string_type keyword = first_token(line, sep);
    
    // heading section
    if (keyword == "heading") {
      
      // print to standard output
      while (ifs.peek() != '*') {
        
        // process line
        if (!getline(ifs, line, lineno))
          continue;
        cout<<line<<endl;
      }
    } // heading section
    
    // node section
    else if (keyword == "node") {
      
      // skip nodes, they will be read by the mesh_io object
      while (ifs.peek() != '*')
        if (!getline(ifs, line, lineno))
          continue;

    } // node section
    
    // node section
    else if (keyword == "element") {
      
      if (!mesh_created) {
        
        // tokenize line
        tokenize(line, tokens);
        
        // get value for key
        key_type eltype = value_for_key("TYPE", tokens);
        
        // get element dimensionality
        Abaqus_element_map abaqus_map;
        
        // create mesh
        mesh = new Mesh(abaqus_map.dimension(eltype));
        
        // object used for reading mesh data
        mesh_io_type mesh_io;

        // read mesh from input file
        mesh_io.read(filename, *mesh);
        
        // set created flag to true
        mesh_created = true;
      }
      
      // skip elements, they will be read by the mesh_io object
      while (ifs.peek() != '*')
        if (!getline(ifs, line, lineno))
          continue;

    } // element section
    
    else {
      
      cout<<"Not implemented keyword "<<keyword<<endl;
      
    }
  }
}


__END_AKANTU__

