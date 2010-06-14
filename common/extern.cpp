/**
 * @file   extern.cpp
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jun 14 14:34:14 2010
 *
 * @brief  initialisation of all global variables
 * to insure the order of creation
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <ostream>

/* -------------------------------------------------------------------------- */
#include "common.hpp"

/* -------------------------------------------------------------------------- */

__BEGIN_MYFEM__

/** \todo write function to get this
 *   values from the environment or a config file
 */

/* -------------------------------------------------------------------------- */
/* error.hpp variables                                                        */
/* -------------------------------------------------------------------------- */
/// standard output for debug messages
std::ostream & _myfem_debug_cout = std::cerr;

/// standard output for normal messages
std::ostream & _myfem_cout = std::cout;

/// debug level
int _debug_level = 6;

/* -------------------------------------------------------------------------- */
/* common.hpp variables                                                       */
/* -------------------------------------------------------------------------- */
int SizeOfType[] = { sizeof(int),
		     sizeof(unsigned int),
		     sizeof(float),
		     sizeof(double)};

const char * TypeCodeName[] = { "int",
				"unsigned short",
				"unsigned int",
				"float",
				"double"};

/* -------------------------------------------------------------------------- */

__END_MYFEM__
