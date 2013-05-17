/**
 * @file   dumper_text.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue May 14 15:27:03 2013
 *
 * @brief  implementation of text dumper
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
#include "dumper_text.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DumperText::DumperText(char separator, bool parallel) : DumperIOHelper() {
  AKANTU_DEBUG_IN();
  
  iohelper::DumperText * dumper_text = new iohelper::DumperText(separator);
  this->dumper = dumper_text;
  
  this->setParallelContext(parallel);

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
