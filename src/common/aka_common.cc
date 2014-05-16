/**
 * @file   aka_common.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Jun 14 19:12:20 2010
 *
 * @brief  Initialization of global variables
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
#include "aka_static_memory.hh"
#include "static_communicator.hh"
#include "aka_random_generator.hh"

#include "parser.hh"
#include "cppargparse.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void initialize(int & argc, char ** & argv) {
  AKANTU_DEBUG_IN();

  initialize("", argc, argv);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void initialize(const std::string & input_file, int & argc, char ** & argv) {
  AKANTU_DEBUG_IN();
  StaticMemory::getStaticMemory();
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator(argc, argv);
  debug::debugger.setParallelContext(comm.whoAmI(), comm.getNbProc());
  debug::initSignalHandler();

  static_argparser.setParallelContext(comm.whoAmI(), comm.getNbProc());
  static_argparser.setExternalExitFunction(debug::Debugger::exit);
  static_argparser.addArgument("--aka_input_file", "Akantu's input file",
			       1, cppargparse::_string, std::string());
  static_argparser.addArgument("--aka_debug_level", "Akantu's overall debug level", 1,
			       cppargparse::_integer, int(dblWarning));

  static_argparser.parse(argc, argv, cppargparse::_remove_parsed);

  std::string infile = static_argparser["aka_input_file"];
  if(infile == "") infile = input_file;

  debug::setDebugLevel(dblError);

  if ("" != infile) {
    static_parser.parse(infile);
  }

  long int seed;
  try {
    seed = static_parser.getParameter("seed", _ppsc_current_scope);
  } catch (debug::Exception & e) {
    seed = time(NULL);
  }

  int dbl_level = static_argparser["aka_debug_level"];
  debug::setDebugLevel(DebugLevel(dbl_level));

  seed *= (comm.whoAmI() + 1);
  Rand48Generator<Real>::seed(seed);
  RandGenerator<Real>::seed(seed);
  AKANTU_DEBUG_INFO("Random seed set to " << seed);

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
void finalize() {
  AKANTU_DEBUG_IN();

  if(StaticMemory::isInstantiated()) delete &(StaticMemory::getStaticMemory());
  if(StaticCommunicator::isInstantiated()) {
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    comm.barrier();
    delete &comm;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
cppargparse::ArgumentParser & getStaticArgumentParser() {
  return static_argparser;
}

/* -------------------------------------------------------------------------- */
const Parser & getStaticParser() {
  return static_parser;
}

/* -------------------------------------------------------------------------- */
const ParserSection & getUserParser() {
  return *(static_parser.getSubSections(_st_user).first);
}

__END_AKANTU__
