/**
 * @file   aka_error.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Sep 06 00:18:58 2010
 *
 * @brief  handling of errors
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
#include "aka_config.hh"
#include "aka_error.hh"
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <csignal>
#include <execinfo.h>
#include <cxxabi.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sys/wait.h>
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
    int status;
    std::string result;
    char * demangled_name;

    if ((demangled_name = abi::__cxa_demangle(symbol, NULL, 0, &status)) != NULL) {
      result = demangled_name;
      free(demangled_name);
    } else {
      result = symbol;
    }

    return result;
    //return symbol;
  }

  /* ------------------------------------------------------------------------ */
  void printBacktrace(__attribute__((unused)) int sig) {
    AKANTU_DEBUG_INFO("Caught  signal " << sig << "!");

    // std::stringstream pidsstr;
    // pidsstr << getpid();
    // char name_buf[512];
    // name_buf[readlink("/proc/self/exe", name_buf, 511)]=0;
    // std::string execname(name_buf);
    // std::cout << "stack trace for " << execname << " pid=" << pidsstr.str() << std::endl;
    // std::string cmd;
    // cmd = "CMDFILE=$(mktemp); echo 'bt' > ${CMDFILE}; gdb --batch " + execname + " " + pidsstr.str() + " < ${CMDFILE};";
    // int retval __attribute__((unused)) = system(("bash -c '" + cmd + "'").c_str());

    const size_t max_depth = 100;
    size_t stack_depth;
    void *stack_addrs[max_depth];
    char **stack_strings;

    size_t i;

    stack_depth = backtrace(stack_addrs, max_depth);
    stack_strings = backtrace_symbols(stack_addrs, stack_depth);

    std::cerr << "BACKTRACE :  " << stack_depth << " stack frames." <<std::endl;
    size_t w = size_t(std::floor(log(double(stack_depth))/std::log(10.))+1);

    /// -1 to remove the call to the printBacktrace function
    for (i = 1; i < stack_depth; i++) {
      std::cerr << "  [" << std::setw(w) << i << "] ";
      std::string bt_line(stack_strings[i]);
      size_t first, second;

      if((first = bt_line.find('(')) != std::string::npos && (second = bt_line.find('+')) != std::string::npos) {
    	std::cerr << bt_line.substr(0,first + 1) << demangle(bt_line.substr(first + 1, second - first - 1).c_str()) <<  bt_line.substr(second) << std::endl;

    	// char name_exe[512];
    	// name_exe[readlink("/proc/self/exe", name_exe, 511)]=0;

    	// std::stringstream syscom;
    	// syscom << "addr2line " << stack_addrs[i] << "-C -s -e " << name_exe;
    	// for (UInt i = 0; i < w + 4; ++i) std::cerr << " ";
    	// std::cout << "-> " << std::flush;
    	// int retval __attribute__((unused)) = system(syscom.str().c_str());
      } else {
    	std::cerr << bt_line << std::endl;
      }
    }

    free(stack_strings);

    std::cerr << "END BACKTRACE" << std::endl;
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
  void Debugger::throwException(const std::string & info) throw(akantu::debug::Exception) {
    AKANTU_DEBUG(akantu::dblWarning, "!!! " << info);
    ::akantu::debug::Exception ex(info, __FILE__, __LINE__ );
    throw ex;
  }

  /* ------------------------------------------------------------------------ */
  void Debugger::exit(int status) {
#ifndef AKANTU_NDEBUG
    int * a = NULL;
    *a = 1;
#endif

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
    UInt pad = std::ceil(std::log10(size));
    sstr << "[R" << std::setfill(' ') << std::right << std::setw(pad)
         << rank << "|S" << size << "] ";
    parallel_context = sstr.str();
  }

  void setDebugLevel(const DebugLevel & level) {
    debugger.setDebugLevel(level);
  }

}

__END_AKANTU__
