/**
 * @file   aka_error.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jun 14 11:43:22 2010
 *
 * @brief  error management and internal exceptions
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
#ifndef __AKANTU_ERROR_HH__
#define __AKANTU_ERROR_HH__

/* -------------------------------------------------------------------------- */
#include <ostream>
#include <sstream>
#include <sys/time.h>

#ifdef AKANTU_USE_MPI
#include <mpi.h>
#endif
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
enum DebugLevel {
  dbl0           = 0,
  dblError       = 0,
  dblAssert      = 0,
  dbl1           = 1,
  dblException   = 1,
  dblCritical    = 1,
  dbl2           = 2,
  dblMajor       = 2,
  dbl3           = 3,
  dblCall        = 3,
  dblSecondary   = 3,
  dblHead        = 3,
  dbl4           = 4,
  dblWarning     = 4,
  dbl5           = 5,
  dblInfo        = 5,
  dbl6           = 6,
  dblIn          = 6,
  dblOut         = 6,
  dbl7           = 7,
  dbl8           = 8,
  dblTrace       = 8,
  dbl9           = 9,
  dblAccessory   = 9,
  dbl10          = 10,
  dblDebug       = 42,
  dbl100         = 100,
  dblDump        = 100,
  dblTest        = 1337
};

/* -------------------------------------------------------------------------- */
namespace debug {
  extern std::ostream & _akantu_cout;

  extern std::ostream *_akantu_debug_cout;

  extern DebugLevel _debug_level;

  extern std::string _parallel_context;

  void setDebugLevel(const DebugLevel & level);
  void setLogFile(const std::string & filename);
  const DebugLevel & getDebugLevel();

  void setParallelContext(int rank, int size);

  void initSignalHandler();
  std::string demangle(const char * symbol);
  void printBacktrace(int sig);


  /* -------------------------------------------------------------------------- */
  /// exception class that can be thrown by akantu
  class Exception : public std::exception {
    /* ------------------------------------------------------------------------ */
    /* Constructors/Destructors                                                 */
    /* ------------------------------------------------------------------------ */
  public:

    //! full constructor
    Exception(std::string info, std::string file, unsigned int line) :
      _info(info), _file(file), _line(line) { }

    //! destructor
    virtual ~Exception() throw() {};

    /* ------------------------------------------------------------------------ */
    /*  Methods                                                                 */
    /* ------------------------------------------------------------------------ */
  public:
    virtual const char* what() const throw() {
      std::stringstream stream;
      stream << "akantu::Exception"
	     << " : " << _info
	     << " ["  << _file << ":" << _line << "]";
      return stream.str().c_str();
    }

    virtual const char* info() const throw() {
      return _info.c_str();
    }

    /* ------------------------------------------------------------------------ */
    /* Class Members                                                            */
    /* ------------------------------------------------------------------------ */
  private:

    /// exception description and additionals
    std::string   _info;

    /// file it is thrown from
    std::string   _file;

    /// ligne it is thrown from
    unsigned int  _line;
  };
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#define AKANTU_LOCATION "(" <<__FILE__ << ":" << __func__ << "():" << __LINE__ << ")"
#ifndef AKANTU_USE_MPI
#define AKANTU_EXIT(status)			\
  do {						\
    int * a = NULL;				\
    *a = 1;					\
    if (status != EXIT_SUCCESS)                 \
      akantu::debug::printBacktrace(15);	\
    exit(status);				\
  } while(0)
#else
#define AKANTU_EXIT(status)			\
  do {						\
    int * a = NULL;				\
    *a = 1;					\
    if (status != EXIT_SUCCESS)                 \
      akantu::debug::printBacktrace(15);	\
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);	\
  } while(0)
#endif

/* -------------------------------------------------------------------------- */
#define AKANTU_EXCEPTION(info)						\
  do {									\
    AKANTU_DEBUG(akantu::dblError, "!!! " << info);			\
    std::stringstream s_info;						\
    s_info << info ;							\
    ::akantu::debug::Exception ex(s_info.str(), __FILE__, __LINE__ );	\
    throw ex;								\
  } while(0)

#define AKANTU_DEBUG_TO_IMPLEMENT()				\
  AKANTU_DEBUG_ERROR(__func__ << " : not implemented yet !")

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_NDEBUG
#define AKANTU_DEBUG_TEST(level)     (false)
#define AKANTU_DEBUG_LEVEL_IS_TEST() (false)
#define AKANTU_DEBUG(level,info)
#define AKANTU_DEBUG_(pref,level,info)
#define AKANTU_DEBUG_IN()
#define AKANTU_DEBUG_OUT()
#define AKANTU_DEBUG_INFO(info)
#define AKANTU_DEBUG_WARNING(info)
#define AKANTU_DEBUG_TRACE(info)
#define AKANTU_DEBUG_ASSERT(test,info)
#define AKANTU_DEBUG_ERROR(info)   AKANTU_EXCEPTION(info)

/* -------------------------------------------------------------------------- */
#else
#define AKANTU_DEBUG(level,info)		\
  AKANTU_DEBUG_("    ",level, info)

#define AKANTU_DEBUG_(pref,level,info)					\
  do {									\
    if(::akantu::debug::_debug_level >= level) {			\
      struct timeval time;						\
      gettimeofday(&time, NULL);					\
      long long  timestamp = time.tv_sec*1e6 + time.tv_usec;/*in us*/	\
      *(::akantu::debug::_akantu_debug_cout)				\
	<< ::akantu::debug::_parallel_context				\
	<< "{" << timestamp << "} "					\
	<< pref << info << " "						\
	<< AKANTU_LOCATION						\
	<< std::endl;							\
    }									\
  } while(0)

#define AKANTU_DEBUG_TEST(level)		\
  (::akantu::debug::_debug_level >= (level))

#define AKANTU_DEBUG_LEVEL_IS_TEST()				\
  (::akantu::debug::_debug_level == (::akantu::dblTest))

#define AKANTU_DEBUG_IN()					\
  AKANTU_DEBUG_("==> ", ::akantu::dblIn,      __func__ << "()")

#define AKANTU_DEBUG_OUT()					\
  AKANTU_DEBUG_("<== ", ::akantu::dblOut,     __func__ << "()")

#define AKANTU_DEBUG_INFO(info)				\
  AKANTU_DEBUG_("--- ", ::akantu::dblInfo,    info)

#define AKANTU_DEBUG_WARNING(info)			\
  AKANTU_DEBUG_("/!\\ ", ::akantu::dblWarning, info)

#define AKANTU_DEBUG_TRACE(info)			\
  AKANTU_DEBUG_(">>> ", ::akantu::dblTrace,   info)

#define AKANTU_DEBUG_ASSERT(test,info)					\
  do {									\
    if (!(test)) {							\
      AKANTU_DEBUG_("!!! ", ::akantu::dblAssert, "assert [" << #test << "] " \
		   << info);						\
      AKANTU_EXIT(EXIT_FAILURE);					\
    }									\
  } while(0)

#define AKANTU_DEBUG_ERROR(info)					\
  do {									\
    AKANTU_DEBUG_("!!! ", ::akantu::dblError, info);			\
    AKANTU_EXIT(EXIT_FAILURE);						\
  } while(0)

#endif // AKANTU_NDEBUG

/* -------------------------------------------------------------------------- */

__END_AKANTU__

#endif /* __AKANTU_ERROR_HH__ */

//  LocalWords:  acessory
