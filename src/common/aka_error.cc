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
#include <cstring>
#include <map>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef AKANTU_USE_MPI
#  include <mpi.h>
#endif

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
  std::string exec(std::string cmd) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "";
    char buffer[1024];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
    result = result.substr(0, result.size()-1);
    pclose(pipe);
    return result;
  }

  /* ------------------------------------------------------------------------ */
  void printBacktrace(__attribute__((unused)) int sig) {
    AKANTU_DEBUG_INFO("Caught  signal " << sig << "!");

#if defined(READLINK_COMMAND) && defined(ADDR2LINE_COMMAND)
    std::string me = "";
    char buf[1024];
    /* The manpage says it won't null terminate.  Let's zero the buffer. */
    memset(buf, 0, sizeof(buf));
    /* Note we use sizeof(buf)-1 since we may need an extra char for NUL. */
    if (readlink("/proc/self/exe", buf, sizeof(buf)-1))
      me = std::string(buf);

    std::ifstream inmaps;
    inmaps.open("/proc/self/maps");
    std::map<std::string, size_t> addr_map;
    std::string line;
    while(inmaps.good()) {
      std::getline(inmaps, line);
      std::stringstream sstr(line);

      size_t first = line.find('-');
      std::stringstream sstra(line.substr(0,first));
      size_t addr; sstra >> std::hex >> addr;

      std::string lib;
      sstr >> lib;
      sstr >> lib;
      sstr >> lib;
      sstr >> lib;
      sstr >> lib;
      sstr >> lib;
      if(lib != "" && addr_map.find(lib) == addr_map.end()) {
        addr_map[lib] = addr;
      }
    }

    if(me != "") addr_map[me] = 0;
#endif

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
      std::cerr << std::dec << "  [" << std::setw(w) << i << "] ";
      std::string bt_line(stack_strings[i]);
      size_t first, second;

      if((first = bt_line.find('(')) != std::string::npos && (second = bt_line.find('+')) != std::string::npos) {
        std::string location = bt_line.substr(0,first);
#if defined(READLINK_COMMAND)
        location = exec(std::string(BOOST_PP_STRINGIZE(READLINK_COMMAND)) + std::string(" -f ") + location);
#endif
        std::string call = demangle(bt_line.substr(first+1, second - first - 1).c_str());
        size_t f = bt_line.find('[');
        size_t s = bt_line.find(']');
        std::string address = bt_line.substr(f+1, s-f-1);
        std::stringstream sstra(address);
        size_t addr; sstra >> std::hex >> addr;

        std::cerr << location << " [" << call << "]";

#if defined(READLINK_COMMAND) && defined(ADDR2LINE_COMMAND)
        std::map<std::string, size_t>::iterator it = addr_map.find(location);
        if(it != addr_map.end()) {
          std::stringstream syscom;
	  syscom << BOOST_PP_STRINGIZE(ADDR2LINE_COMMAND) << " 0x" << std::hex << (addr - it->second) << " -i -e " << location;
          std::string line = exec(syscom.str());
          std::cerr << " (" << line << ")" <<std::endl;
        } else {
#endif
          std::cerr << " (0x" << std::hex << addr << ")" <<std::endl;
#if defined(READLINK_COMMAND) && defined(ADDR2LINE_COMMAND)
        }
#endif
      } else {
    	std::cerr << bt_line << std::endl;
      }
    }

    free(stack_strings);

    std::cerr << "END BACKTRACE" << std::endl;
    if(sig != 15) debugger.exit(EXIT_FAILURE);
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
    AKANTU_DEBUG_( "###", akantu::dblWarning, info);
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
    sstr << "<" << getpid() << ">[R" << std::setfill(' ') << std::right << std::setw(pad)
         << rank << "|S" << size << "] ";
    parallel_context = sstr.str();
  }

  void setDebugLevel(const DebugLevel & level) {
    debugger.setDebugLevel(level);
  }

}

__END_AKANTU__
