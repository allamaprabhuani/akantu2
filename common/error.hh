/**
 * @file   error.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jun 14 11:43:22 2010
 *
 * @brief  error management and internal exceptions
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_ERROR_HH__
#define __AKANTU_ERROR_HH__

/* -------------------------------------------------------------------------- */
#include <ostream>
#include <sstream>

/* -------------------------------------------------------------------------- */
#include "common.hh"

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
  dblDump        = 100
};
/* -------------------------------------------------------------------------- */

extern std::ostream & _akantu_cout;

extern std::ostream & _akantu_debug_cout;

extern DebugLevel _debug_level;

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

/* -------------------------------------------------------------------------- */

#define AKANTU_LOCATION "(" <<__FILE__ << ":" << __func__ << "():" << __LINE__ << ")"
#ifndef AKANTU_MPI
#define AKANTU_EXIT(status)			\
  do {						\
    exit(status);				\
  } while(0)
#else
#define AKANTU_EXIT(status)			\
  do {						\
    MPI_Finalize();				\
    exit(status);				\
  } while(0)
#endif
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_NDEBUG
#define AKANTU_DEBUG_TEST(level)   (0)
#define AKANTU_DEBUG(level,info)
#define AKANTU_DEBUG_IN()
#define AKANTU_DEBUG_OUT()
#define AKANTU_DEBUG_INFO(info)
#define AKANTU_DEBUG_WARNING(info)
#define AKANTU_DEBUG_TRACE(info)
#define AKANTU_DEBUG_ASSERT(test,info)
#define AKANTU_DEBUG_ERROR(info)					\
  do {									\
    AKANTU_DEBUG(akantu::dblError, "!!! " << info);			\
    std::stringstream s_info;						\
    s_info << info ;							\
    Exception ex(s_info.str(), __FILE__, __LINE__ );			\
    throw ex;								\
  } while(0)

/* -------------------------------------------------------------------------- */
#else
#define AKANTU_DEBUG(level,info)					\
  ((akantu::_debug_level >= level) &&					\
   (_akantu_debug_cout << info << " " << AKANTU_LOCATION << std::endl))

#define AKANTU_DEBUG_TEST(level)		\
  (akantu::_debug_level >= (level))

#define AKANTU_DEBUG_IN()					\
  AKANTU_DEBUG(akantu::dblIn     , "==> " << __func__ << "()")

#define AKANTU_DEBUG_OUT()					\
  AKANTU_DEBUG(akantu::dblOut    , "<== " << __func__ << "()")

#define AKANTU_DEBUG_INFO(info)				\
  AKANTU_DEBUG(akantu::dblInfo   , "--- " << info)

#define AKANTU_DEBUG_WARNING(info)			\
  AKANTU_DEBUG(akantu::dblWarning, "??? " << info)

#define AKANTU_DEBUG_TRACE(info)			\
  AKANTU_DEBUG(akantu::dblTrace  , ">>> " << info)

#define AKANTU_DEBUG_ASSERT(test,info)					\
  do {									\
    if (!(test)) {							\
      AKANTU_DEBUG(dblAssert, "assert [" << #test << "] "		\
		   << "!!! " << info);					\
      AKANTU_EXIT(EXIT_FAILURE);					\
    }									\
  } while(0)

#define AKANTU_DEBUG_ERROR(info)					\
  do {									\
    AKANTU_DEBUG(akantu::dblError, "!!! " << info);			\
    AKANTU_EXIT(EXIT_FAILURE);						\
  } while(0)

#endif

/* -------------------------------------------------------------------------- */

__END_AKANTU__

#endif /* __AKANTU_ERROR_HH__ */

//  LocalWords:  acessory
