/**
 * @file   file_manager.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Dec 06 2012
 *
 * @brief  file manager header
 *
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under
 * the terms  of the  GNU Lesser  General Public  License as  published by  the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for
 * more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef IOHELPER_FILE_MANAGER_H_
#define IOHELPER_FILE_MANAGER_H_
/* -------------------------------------------------------------------------- */
#include "iohelper_common.hh"
#include <fstream>
#include <iostream>
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <typeinfo>
/* -------------------------------------------------------------------------- */

namespace iohelper {

template <class charT, class Traits = std::char_traits<charT>>
class GZfstream : public std::basic_fstream<charT, Traits> {

public:
  GZfstream();
  GZfstream(const std::string & fname,
            std::fstream::openmode mode = std::fstream::out,
            bool compr = false);

  ///! opening methods
  inline void open(const std::string & name,
                   std::fstream::openmode mode = std::fstream::out,
                   bool compr = false);

private:
};

template <class charT, class Traits>
inline GZfstream<charT, Traits>::GZfstream()
    : std::basic_fstream<charT, Traits>() {}

template <class charT, class Traits>
inline GZfstream<charT, Traits>::GZfstream(const std::string & fname,
                                           std::fstream::openmode mode,
                                           bool /*unused*/)
    : std::basic_fstream<charT, Traits>(fname.c_str(), mode) {}

template <class charT, class Traits>
inline void GZfstream<charT, Traits>::open(const std::string & fname,
                                           std::fstream::openmode mode,
                                           bool /*unused*/) {
  std::basic_fstream<charT, Traits>::open(fname.c_str(), mode);
}

using File = GZfstream<char>;

} // namespace iohelper

#endif /* IOHELPER_FILE_MANAGER_H_ */
