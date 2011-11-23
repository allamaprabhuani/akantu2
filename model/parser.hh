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

__BEGIN_AKANTU__

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

  /// open a file to parse
  void open(const std::string & filename);

  /// read the file and return the next material type
  std::string getNextSection(const std::string & obj_type);

  /// read properties and instanciate a given material object
  template <typename M>
  void readSection(M & model);
  template <typename Obj>
  Obj * readSection(Model & model, std::string & obj_name);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  __aka_inline__ void my_getline();

  /// number of the current line
  UInt current_line;

  /// material file
  std::ifstream infile;

  /// current read line
  std::string line;
};


#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "parser_inline_impl.cc"
#endif

__END_AKANTU__

#endif
