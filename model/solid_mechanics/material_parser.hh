/**
 * @file   material_parser.hh
 * @author Guillaume ANCIAUX <anciaux@lsmscluster1.epfl.ch>
 * @date   Thu Nov 25 11:43:48 2010
 *
 * @brief  
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_PARSER_HH__
#define __AKANTU_MATERIAL_PARSER_HH__

__BEGIN_AKANTU__

class MaterialParser {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  MaterialParser(){};
  virtual ~MaterialParser(){infile.close();};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// open a file to parse
  void open(const std::string & filename){
    infile.open(filename.c_str());
    
    current_line = 0;
    
    if(!infile.good()) {
      AKANTU_DEBUG_ERROR("Cannot open file " << filename);
    }
  };

  /// read the file and return the next material type
  std::string getNextMaterialType(){
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
  };

  template <typename M> 
  Material * readMaterialObject(SolidMechanicsModel & model,MaterialID & mat_id){
    std::string keyword;
    std::string value;
    
    /// instanciate the material object
    Material * mat = new M(model, mat_id);
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
      } catch (Exception ex) {
  	AKANTU_DEBUG_ERROR("Malformed material file : error in setParam \""
  			   << ex.info() << "\" at line " << current_line);
      }
      
      my_getline();
    }
    return mat;
  };

  
  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:


  inline void my_getline() {
    std::getline(infile, line); //read the line
    if (!(infile.flags() & (std::ios::failbit | std::ios::eofbit)))
      ++current_line;
    size_t pos = line.find("#"); //remove the comment
    line = line.substr(0, pos);
    trim(line); // remove unnecessary spaces
  }

  UInt current_line;
  std::ifstream infile;
  std::string line;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_parser_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const MaterialParser & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif
