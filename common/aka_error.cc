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
#include "aka_error.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

namespace debug {
  void setDebugLevel(const DebugLevel & level) {
    _debug_level = level;
  }

  void setParallelContext(int rank, int size) {
    std::stringstream sstr;
    sstr << "[" << std::setfill(' ') << std::right << std::setw(3) << (rank + 1) << "/" << size << "] ";
    _parallel_context = sstr.str();
  }
};

__END_AKANTU__
