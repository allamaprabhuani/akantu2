/**
 * @file   dumper.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Mon Oct 02 2017
 *
 * @brief Dumper interface
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
#include "aka_array.hh"
#include "element_type_map.hh"
#include "support.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <memory>
#include <string>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DUMPER_HH
#define AKANTU_DUMPER_HH

namespace akantu {

// TODO: rename class hierarchy since they also handle read
class Dumper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit Dumper(dumper::SupportBase & support);
  virtual ~Dumper();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// request dump: this calls the internal dumper implementation
  void dump();

  /// request dump: this calls the internal dumper implementation
  void read();

  /// set the directory where to generate the dumped files
  virtual void setDirectory(const std::string & directory);
  /// set the base name (needed by most IOHelper dumpers)
  virtual void setBaseName(const std::string & basename);

protected:
  virtual void dumpInternal() = 0;

  virtual void readInternal() { AKANTU_TO_IMPLEMENT(); }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Set the absolute time of the next dump
  void setTime(Real time) {
    this->time = time;
    this->time_activated = true;
  }

  /// Set the amount to increase the time at each dump
  void setTimeStep(Real time_step) {
    this->time_step = time_step;
    this->time_activated = true;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  dumper::SupportBase & support;

  /// directory name
  std::string directory{"./"};

  /// filename prefix
  std::string filename{"dump"};

  /// dump counter
  UInt count{0};

  /// is time tracking activated in the dumper
  bool time_activated{true};

  Real time{0.};
  Real time_step{0.};

  Int prank{0};
  Int psize{1};
};

} // namespace akantu

#endif /* AKANTU_DUMPER_HH */
