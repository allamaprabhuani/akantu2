/**
 * @file   material_parser_inline_impl.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Nov 26 08:32:38 2010
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
inline std::string Parser::getNextSection(const std::string & section_type){
  while(infile.good()) {
    my_getline();

    // if empty line continue
    if(line.empty()) continue;

    std::stringstream sstr(line);
    std::string keyword;
    std::string value;

    sstr >> keyword;
    to_lower(keyword);
    /// if found a material deccription then stop
    /// and prepare the things for further reading
    if(keyword == section_type) {
      std::string type; sstr >> type;
      to_lower(type);
      std::string obracket; sstr >> obracket;
      if(obracket != "[")
	AKANTU_DEBUG_ERROR("Malformed config file : missing [ at line " << current_line);
      return type;
    }
  }
  return "";
}

/* -------------------------------------------------------------------------- */
inline void Parser::my_getline() {
  std::getline(infile, line); //read the line
  if (!(infile.flags() & (std::ios::failbit | std::ios::eofbit)))
    ++current_line;
  size_t pos = line.find("#"); //remove the comment
  line = line.substr(0, pos);
  trim(line); // remove unnecessary spaces
  std::cout << line << std::endl;
}

/* -------------------------------------------------------------------------- */
inline void Parser::open(const std::string & filename){
  infile.open(filename.c_str());
  current_line = 0;

  if(!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }
}
