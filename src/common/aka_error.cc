/**
 * @file   aka_error.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Sun Sep  5 21:03:37 2010
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
#include <csignal>
//#include <cerrno>
#include <execinfo.h>
//#include <cxxabi.h>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_error.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

namespace debug {
  /* ------------------------------------------------------------------------ */
  void initSignalHandler() {
    struct sigaction action;

    action.sa_handler = &printBacktrace;
    sigemptyset(&(action.sa_mask));
    action.sa_flags = SA_RESETHAND;

    sigaction(SIGSEGV, &action, NULL);
  }

  /* ------------------------------------------------------------------------ */
  std::string demangle(const char* symbol) {
    // int status;
    // std::string result;
    // char * demangled_name;

    // if ((demangled_name = abi::__cxa_demangle(symbol, NULL, 0, &status)) != NULL) {
    //   result = demangled_name;
    //   free(demangled_name);
    // } else {
    //   result = symbol;
    // }

    // return result;
    return symbol;
  }


  /* ------------------------------------------------------------------------ */
  void printBacktrace(int sig) {
    AKANTU_DEBUG_INFO("Caught  signal " << sig << "!");

    void *array[10];
    size_t size;
    char **strings;
    size_t i;

    size = backtrace (array, 10);
    strings = backtrace_symbols (array, size);

    std::cerr << "BACKTRACE :  " << size - 1 << " stack frames." <<std::endl;
    /// -1 to remove the call to the printBacktrace function
    for (i = 1; i < size; i++)
      std::cerr << "   " << demangle(strings[i]) << std::endl;;

    free (strings);

    std::cerr << "END BACKTRACE" << std::endl;

    // int * segfault = NULL;
    // *segfault = 0;
  }

  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */

  Debugger::Debugger() {
    cout = &std::cerr;
    level = dblInfo;
    parallel_context = "";
    file_open = false;
  }

  /* ------------------------------------------------------------------------ */
  Debugger::~Debugger() {
    if(file_open) {
      dynamic_cast<std::ofstream *>(cout)->close();
      delete cout;
    }
  }

  /* ------------------------------------------------------------------------ */
  void Debugger::throwException(const std::string & info) {
    AKANTU_DEBUG(akantu::dblWarning, "!!! " << info);
    ::akantu::debug::Exception ex(info, __FILE__, __LINE__ );
    throw ex;
  }

  /* ------------------------------------------------------------------------ */
  void Debugger::exit(int status) {
    int * a = NULL;
    *a = 1;
    if (status != EXIT_SUCCESS)
      akantu::debug::printBacktrace(15);
#ifdef AKANTU_USE_MPI
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
#endif
    exit(status); // not  called when compiled  with MPI  due to  MPI_Abort, but
                  // MPI_Abort does not have the noreturn attribute
  }

  /* ------------------------------------------------------------------------ */
  void Debugger::setDebugLevel(const DebugLevel & level) {
    this->level = level;
  }
  /* ------------------------------------------------------------------------ */
  const DebugLevel & Debugger::getDebugLevel() const {
    return level;
  }
  /* ------------------------------------------------------------------------ */
  void Debugger::setLogFile(const std::string & filename) {
    if(file_open) {
      dynamic_cast<std::ofstream *>(cout)->close();
      delete cout;
    }
    std::ofstream * fileout = new std::ofstream(filename.c_str());
    file_open = true;
    cout = fileout;
  }

  std::ostream & Debugger::getOutputStream() {
    return *cout;
  }

  /* ------------------------------------------------------------------------ */
  void Debugger::setParallelContext(int rank, int size) {
    std::stringstream sstr;
    sstr << "[" << std::setfill(' ') << std::right << std::setw(3)
         << (rank + 1) << "/" << size << "] ";
    parallel_context = sstr.str();
  }

  void setDebugLevel(const DebugLevel & level) {
    debugger.setDebugLevel(level);
  }

}

__END_AKANTU__
