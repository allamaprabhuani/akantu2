/**
 * @file   parser_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Nov 24 09:36:33 2011
 *
 * @brief  Template part of the parser
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

template <typename M>
inline void Parser::readSection(M & model) {
  my_getline();

  while(line[0] != ']') {
    if(line != "") {
      size_t pos = line.find("=");
      if(pos == std::string::npos)
	AKANTU_DEBUG_ERROR("Malformed material file : line must be \"key = value\" at line"
			   << current_line);

      std::string keyword = trim(line.substr(0, pos));
      std::string value   = trim(line.substr(pos + 1));

      if(!model.setParam(keyword, value)) {
	AKANTU_DEBUG_ERROR("Malformed file : error in setParam at line "
			   << current_line <<"."
			   << " Parameter (" << keyword << ") is not recognized!");
      }
    }
    my_getline();
  }
}

/* -------------------------------------------------------------------------- */
inline void Parser::readSection(const std::string & obj_name,
				Parsable & obj){
  /// read the material properties
  my_getline();

  while(line[0] != ']') {
    if(line != "") {
      size_t pos = line.find("=");
      if(pos == std::string::npos)
	AKANTU_DEBUG_ERROR("Malformed material file : line must be \"key = value\" at line"
			   << current_line);

      std::string keyword = trim(line.substr(0, pos));
      std::string value   = trim(line.substr(pos + 1));

      if(!obj.parseParam(keyword, value, obj_name)) {
	AKANTU_DEBUG_WARNING("Malformed material file : error in setParam at line "
			     << current_line <<"."
			     << " Parameter (" << keyword << ") is not recognized by " << obj_name << "!");
      }
    }
    my_getline();
  }
}

