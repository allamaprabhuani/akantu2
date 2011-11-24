/**
 * @file   parser_tmpl.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Nov 24 08:44:22 2011
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

/* -------------------------------------------------------------------------- */
template <typename M>
inline void Parser::readSection(M & model){
  std::string keyword;
  std::string value;

  my_getline();

  while(line[0] != ']') {
    if (line.empty()){
      my_getline();
      continue;
    }
    size_t pos = line.find("=");
    if(pos == std::string::npos)
      AKANTU_DEBUG_ERROR("Malformed config file : line must be \"key = value\" at line"
			 << current_line);

    keyword = line.substr(0, pos);  trim(keyword);
    value   = line.substr(pos + 1); trim(value);

    try {
      model.setParam(keyword, value);
    } catch (debug::Exception ex) {
      AKANTU_DEBUG_ERROR("Malformed config file : error in setParam \""
			 << ex.info() << "\" at line " << current_line);
    }

    my_getline();
  }
}

/* -------------------------------------------------------------------------- */
template <typename Obj>
inline Obj * Parser::readSection(Model & model,std::string & obj_name){
  std::string keyword;
  std::string value;

  /// instanciate the material object
  Obj * obj = new Obj(model, obj_name);
  /// read the material properties
  my_getline();

  while(line[0] != ']') {
    size_t pos = line.find("=");
    if(pos == std::string::npos)
      AKANTU_DEBUG_ERROR("Malformed config file : line must be \"key = value\" at line"
			 << current_line);

    keyword = line.substr(0, pos);  trim(keyword);
    value   = line.substr(pos + 1); trim(value);

    try {
      obj->setParam(keyword, value, obj_name);
    } catch (debug::Exception ex) {
      AKANTU_DEBUG_ERROR("Malformed config file : error in setParam \""
			 << ex.info() << "\" at line " << current_line);
    }

    my_getline();
  }
  return obj;
}
