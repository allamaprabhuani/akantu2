/**
 * @file   material_parser.hh
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Thu Nov 25 11:43:48 2010
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

#ifndef __AKANTU_MATERIAL_PARSER_HH__
#define __AKANTU_MATERIAL_PARSER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */
#include <fstream>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Parsable {
public:
  virtual ~Parsable() {};
  virtual bool parseParam(const std::string & key, const std::string & value,
			  const ID & id) = 0;
};

class Parser {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Parser() : current_line(0) {};
  virtual ~Parser(){ infile.close(); };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// open/close a file to parse
  void open(const std::string & filename);
  void close();

  /// read the file and return the next material type
  std::string getNextSection(const std::string & obj_type, std::string & optional_param);

  /// read properties and instanciate a given material object
  template <typename M>
  void readSection(M & model);

  /// read the properties in a section
  void readSection(const std::string & obj_name, Parsable & obj);

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

  /// name of file parsed
  std::string filename;
};

#include "parser_tmpl.hh"

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "parser_inline_impl.cc"
#endif

__END_AKANTU__

#endif
