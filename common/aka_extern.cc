/**
 * @file   aka_extern.cpp
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jun 14 14:34:14 2010
 *
 * @brief  initialisation of all global variables
 * to insure the order of creation
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <ostream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/** \todo write function to get this
 *   values from the environment or a config file
 */

/* -------------------------------------------------------------------------- */
/* error.hpp variables                                                        */
/* -------------------------------------------------------------------------- */
/// standard output for debug messages
std::ostream & _akantu_debug_cout = std::cerr;

/// standard output for normal messages
std::ostream & _akantu_cout = std::cout;

/// debug level
DebugLevel _debug_level = (DebugLevel) 5;

/* -------------------------------------------------------------------------- */

__END_AKANTU__
