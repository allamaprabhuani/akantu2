/**
 * @file   aka_abaqus_parser.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Fri Oct  7 10:47:00 2011
 *
 * @brief  Objects and functions needed for parsing Abaqus input files
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

#ifndef __AKANTU_ABAQUS_PARSER_HH__
#define __AKANTU_ABAQUS_PARSER_HH__

#include <iostream>

// other external library header files
#include<boost/tokenizer.hpp>

__BEGIN_AKANTU__

using std::cout;
using std::endl;

//! Check for empty string
template <class string_type>
bool empty_string(const string_type& s) {
  
  // loop over string characters
  for (int i=0; i<s.size(); ++i)
    // check if character is a space character: ' ' space, '\t' horizontal tab,
    // '\n' newline, '\v' vertical tab, '\f' feed, '\r' carriage return
    if (!isspace(s[i]))
      return false;
  return true;
}


template <class input_stream_type, class string_type, typename counter_type>
bool getline(input_stream_type& is, string_type& str, counter_type& n) {
  
  // get line from file input stream
  std::getline(is, str);
  
  // increment line counter
  ++n;
  
  // get position of comment sequence
  size_t pos = str.find("**");
  if (pos == 0)
    return false;
  else if (pos != std::string::npos)
    str = str.substr(pos, str.size());
  // else line is not a comment
  
  // skip line if empty
  if (empty_string(str)) {
    cout<<"*** WARNING *** Found empty line at location "<<n<<endl;
    return false;
  }
  return true;
}


template <class container_type>
void tokenize(const typename container_type::value_type& str,
              container_type& tokens) {
  
  typedef typename container_type::value_type string_type;
  typedef typename string_type::value_type value_type;
  
  boost::char_separator<value_type> sep(" \t*=,\r\n");
  tokenize(str, tokens, sep);
}


template <class container_type, class separator_list>
void tokenize(const typename container_type::value_type& str,
              container_type& tokens,
              separator_list& sep) {
  
  typedef typename container_type::value_type string_type;
  typedef typename string_type::value_type value_type;
  typedef typename string_type::traits_type traits_type;
  typedef typename string_type::const_iterator const_iterator;
  typedef boost::tokenizer< separator_list, const_iterator, string_type> tokenizer_type;
  
  // tokenize string
  tokenizer_type t(str, sep);
  
  // add tokens to container
  tokens.assign(t.begin(), t.end());
}

template <class string_type, class separator_list>
string_type first_token(const string_type& str, separator_list& sep) {
  
  typedef typename string_type::value_type value_type;
  typedef typename string_type::traits_type traits_type;
  typedef typename string_type::const_iterator const_iterator;
  typedef boost::tokenizer< separator_list, const_iterator, string_type> tokenizer_type;
  
  // tokenize string
  tokenizer_type t(str, sep);
  
  return *t.begin();
}


template <class container_type>
const typename container_type::value_type& value_for_key(const char* key, const container_type& tokens) {
  
  for (typename container_type::const_iterator it = tokens.begin(); it != tokens.end(); ++it)
    if (*it == key) {
      // advance and return dereferenced iterator
      assert (++it != tokens.end());
      return *it;
    }
}


template <typename T, class string_type, typename counter_type>
T stream_cast(const string_type& s, counter_type n) {
  
  typedef typename string_type::value_type char_type;
  typedef typename string_type::traits_type traits_type;
  
  typename std::basic_istringstream<char_type, traits_type> iss(s);
  
  T x = T();
  char c;
  
  iss >> x;
//  cout<<" iss bad -> "<<iss.bad()<<endl;


  if (!(iss) || iss.get(c)) {
    cout<<"*** WARNING *** Bad conversion to type "<<typeid(T).name()<<" in line "<<n<<endl;
    exit(EXIT_FAILURE);
  }
  return x;
}


template <typename T, typename counter_type>
T stream_cast(const char* cstr, counter_type n) {
  return stream_cast<T>(std::string(cstr), n);
}
  

// \todo remove this function after finding out why the badbit is set in 
// ci_string
template <typename T, typename counter_type>
T stream_cast(const ci_string& str, counter_type n) {
  return stream_cast<T>(str.data(), n);
}




class Abaqus_element_map {
  
  typedef ci_string key_type;
  typedef std::map<key_type, ElementType> map_type;
  typedef map_type::const_iterator map_iterator;
  
  map_type map;
  
public:
  
  Abaqus_element_map() {
    
    // 2D element types
    
    map["CPE3"] = _triangle_3;
    map["CPS3"] = _triangle_3;
    map["DC2D3"] = _triangle_3;
    
    map["CPE6"] = _triangle_6;
    map["CPS6"] = _triangle_6;
    map["DC2D6"] = _triangle_6;
    
    map["CPE4"] = _quadrangle_4;
    map["CPS4"] = _quadrangle_4;
    map["DC2D4"] = _quadrangle_4;
    
    map["CPE8"] = _quadrangle_8;
    map["CPS8"] = _quadrangle_8;
    map["DC2D8"] = _quadrangle_8;
    
    // 3D element types
    
    map["C3D4"] = _tetrahedron_4;
    map["DC3D4"] = _tetrahedron_4;
    
    map["C3D8"] = _hexahedron_8;
    map["DC3D8"] = _hexahedron_8;
    
    map["C3D10"] = _tetrahedron_10;
    map["DC3D10"] = _tetrahedron_10;
  }
  
  ElementType operator[](const key_type& id) {
    
    map_iterator it = map.find(id);
    if (it == map.end()) {
      cout<<"*** ERROR *** Abaqus element with ID "<<id<<" not found in map"<<endl;
      cout<<"*** ABORT ***"<<endl;
      exit(EXIT_FAILURE);
    }
    return it->second;
  }
  
  UInt dimension(const key_type& id) {
    
    map_iterator it = map.find(id);
    if (it == map.end()) {
      cout<<"*** ERROR *** Abaqus element with ID "<<id<<" not found in map"<<endl;
      cout<<"*** ABORT ***"<<endl;
      exit(EXIT_FAILURE);
    }
    return Mesh::getSpatialDimension(it->second);
  }
  
};


__END_AKANTU__

#endif /* __AKANTU_ABAQUS_PARSER_HH__ */
