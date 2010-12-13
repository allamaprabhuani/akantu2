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
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
template <typename Mat>
inline Material * MaterialParser::readMaterialObject(SolidMechanicsModel & model,MaterialID & mat_id){
  std::string keyword;
  std::string value;

  /// instanciate the material object
  Material * mat = new Mat(model, mat_id);
  /// read the material properties
  my_getline();

  while(line[0] != ']') {
    size_t pos = line.find("=");
    if(pos == std::string::npos)
      AKANTU_DEBUG_ERROR("Malformed material file : line must be \"key = value\" at line"
			 << current_line);

    keyword = line.substr(0, pos);  trim(keyword);
    value   = line.substr(pos + 1); trim(value);

    try {
      mat->setParam(keyword, value, mat_id);
    } catch (debug::Exception ex) {
      AKANTU_DEBUG_ERROR("Malformed material file : error in setParam \""
			 << ex.info() << "\" at line " << current_line);
    }

    my_getline();
  }
  return mat;
}

/* -------------------------------------------------------------------------- */
inline std::string MaterialParser::getNextMaterialType(){
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
    if(keyword == "material") {
      std::string type; sstr >> type;
      to_lower(type);
      std::string obracket; sstr >> obracket;
      if(obracket != "[")
	AKANTU_DEBUG_ERROR("Malformed material file : missing [ at line " << current_line);
      return type;
    }
  }
  return "";
}

/* -------------------------------------------------------------------------- */
inline void MaterialParser::my_getline() {
  std::getline(infile, line); //read the line
  if (!(infile.flags() & (std::ios::failbit | std::ios::eofbit)))
    ++current_line;
  size_t pos = line.find("#"); //remove the comment
  line = line.substr(0, pos);
  trim(line); // remove unnecessary spaces
}

/* -------------------------------------------------------------------------- */
inline void MaterialParser::open(const std::string & filename){
  infile.open(filename.c_str());
  current_line = 0;

  if(!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }
}
