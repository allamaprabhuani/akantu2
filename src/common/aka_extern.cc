/**
 * @file   aka_extern.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Jun 14 19:12:20 2010
 *
 * @brief  initialisation of all global variables
 * to insure the order of creation
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
#include "aka_common.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/** \todo write function to get this
 *   values from the environment or a config file
 */

/* -------------------------------------------------------------------------- */
/* error.hpp variables                                                        */
/* -------------------------------------------------------------------------- */
namespace debug {
  /// standard output for debug messages
  std::ostream *_akantu_debug_cout = &std::cerr;

  /// standard output for normal messages
  std::ostream & _akantu_cout = std::cout;

  /// debug level
  DebugLevel _debug_level = dblInfo;

  /// parallel context used in debug messages
  std::string _parallel_context = "";

  Debugger debugger;

}

Real Math::tolerance = 1e-8;

/* -------------------------------------------------------------------------- */

__END_AKANTU__
