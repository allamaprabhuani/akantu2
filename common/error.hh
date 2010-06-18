/**
 * @file   error.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jun 14 11:43:22 2010
 *
 * @brief  error management and internal exceptions
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __MYFEM_ERROR_HH__
#define __MYFEM_ERROR_HH__

/* -------------------------------------------------------------------------- */
#include <ostream>
#include <sstream>

/* -------------------------------------------------------------------------- */
#include "common.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_MYFEM__

/* -------------------------------------------------------------------------- */
extern std::ostream & _myfem_cout;

extern std::ostream & _myfem_debug_cout;

extern int _debug_level;

/* -------------------------------------------------------------------------- */
/// exception class that can be thrown by myfem
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
    stream << "myfem::Exception"
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
enum DebugLevel {
  dbl0           = 0,
  dblError       = 0,
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
  dbl10          = 10,
  dblDump        = 100
};


#define MYFEM_LOCATION "(" <<__FILE__ << ":" << __LINE__ << ") "
#ifndef MYFEM_MPI
#define MYFEM_EXIT(status) exit(status)
#else
#define MYFEM_EXIT(status)			\
  do {						\
    MPI_Finalize();				\
    exit(status);				\
  } while(0)
#endif
/* -------------------------------------------------------------------------- */
#ifdef MYFEM_NDEBUG
#define MYFEM_DEBUG_TEST(level)   (0)
#define MYFEM_DEBUG(level,info)
#define MYFEM_DEBUG_IN()
#define MYFEM_DEBUG_OUT()
#define MYFEM_DEBUG_INFO(info)
#define MYFEM_DEBUG_WARNING(info)
#define MYFEM_DEBUG_TRACE(info)
#define MYFEM_DEBUG_ASSERT(test,info)
/* -------------------------------------------------------------------------- */
#else
#define MYFEM_DEBUG(level,info)						\
  ((myfem::_debug_level >= level) &&					\
   (_myfem_debug_cout << info << " " << MYFEM_LOCATION << std::endl))

#define MYFEM_DEBUG_TEST(level)			\
  (myfem::_debug_level >= (level))

#define MYFEM_DEBUG_IN()					\
  MYFEM_DEBUG(myfem::dblIn     , "==> " << __func__ << "()")

#define MYFEM_DEBUG_OUT()					\
  MYFEM_DEBUG(myfem::dblOut    , "<== " << __func__ << "()")

#define MYFEM_DEBUG_INFO(info)				\
  MYFEM_DEBUG(myfem::dblInfo   , "--- " << info)

#define MYFEM_DEBUG_WARNING(info)			\
  MYFEM_DEBUG(myfem::dblWarning, "??? " << info)

#define MYFEM_DEBUG_TRACE(info)				\
  MYFEM_DEBUG(myfem::dblTrace  , ">>> " << info)

#define MYFEM_DEBUG_ASSERT(test,info)					\
  do {									\
    if (!(test)) {							\
      _myfem_debug_cout << "(" <<__FILE__ << ":" << __LINE__ << ") "	\
			<< "assert [" << #test << "] "			\
			<< "!!! " << info				\
			<< std::endl;					\
      MYFEM_EXIT(EXIT_FAILURE);						\
    }									\
  } while(0)
#endif

/* -------------------------------------------------------------------------- */
#define MYFEM_DEBUG_ERROR(info)						\
  do {									\
    MYFEM_DEBUG(myfem::dblError, "!!! " << info);			\
    std::stringstream s_info;						\
    s_info << info ;							\
    Exception ex(s_info.str(), __FILE__, __LINE__ );			\
    throw ex;								\
  } while(0)

/* -------------------------------------------------------------------------- */

__END_MYFEM__

#endif /* __MYFEM_ERROR_HH__ */
