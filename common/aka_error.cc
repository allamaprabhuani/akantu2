/**
 * @file   aka_error.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Sun Sep  5 21:03:37 2010
 *
 * @brief
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <csignal>
#include <cerrno>
#include <execinfo.h>
#include <cxxabi.h>

/* -------------------------------------------------------------------------- */
#include "aka_error.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

namespace debug {
  /* ------------------------------------------------------------------------ */
  void setDebugLevel(const DebugLevel & level) {
    _debug_level = level;
  }

  /* ------------------------------------------------------------------------ */
  void setParallelContext(int rank, int size) {
    std::stringstream sstr;
    sstr << "[" << std::setfill(' ') << std::right << std::setw(3) << (rank + 1) << "/" << size << "] ";
    _parallel_context = sstr.str();
  }

  /* ------------------------------------------------------------------------ */
  void initSignalHandler() {
    struct sigaction action;

    action.sa_handler = &printBacktrace;
    sigemptyset(&(action.sa_mask));
    action.sa_flags = SA_RESETHAND;

    sigaction(SIGSEGV, &action, NULL);
  }

  /* ------------------------------------------------------------------------ */
  static std::string demangle(const char* symbol) {
    std::string trace(symbol);

    std::string::size_type begin = trace.find_first_of('(') + 1;
    std::string::size_type end   = trace.find_last_of('+');

    if (begin != std::string::npos && end != std::string::npos) {
      std::string sub_trace = trace.substr(begin, end - begin);

      int status;
      size_t size;
      std::string result;
      char * demangled_name;

      if ((demangled_name = abi::__cxa_demangle(sub_trace.c_str(), NULL, &size, &status)) != NULL) {
	result = demangled_name;
	free(demangled_name);
      }

      std::stringstream result_sstr;
      result_sstr << trace.substr(0, begin) <<  result << trace.substr(end);

      return result_sstr.str();
    }

    return trace;
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

    int * segfault = NULL;
    *segfault = 0;
  }

};

__END_AKANTU__
