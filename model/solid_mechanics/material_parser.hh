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

#ifndef __AKANTU_MATERIAL_PARSER_HH__
#define __AKANTU_MATERIAL_PARSER_HH__

__BEGIN_AKANTU__

class MaterialParser {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialParser() : current_line(0) {};
  virtual ~MaterialParser(){ infile.close(); };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// open a file to parse
  void open(const std::string & filename);

  /// read the file and return the next material type
  std::string getNextMaterialType();

  /// read properties and instanciate a given material object
  template <typename Mat>
  Material * readMaterialObject(SolidMechanicsModel & model, MaterialID & mat_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  inline void my_getline();

  /// number of the current line
  UInt current_line;

  /// material file
  std::ifstream infile;

  /// current read line
  std::string line;
};


#include "material_parser_inline_impl.cc"

__END_AKANTU__

#endif
