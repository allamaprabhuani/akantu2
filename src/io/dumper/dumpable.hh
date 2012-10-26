/**
 * @file   dumpable.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Oct 24 22:37:42 2012
 *
 * @brief  Interface for object who wants to dump themselves
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

#ifndef __AKANTU_DUMPABLE_HH__
#define __AKANTU_DUMPABLE_HH__

#ifdef AKANTU_USE_IOHELPER
#include "dumper_paraview.hh"
__BEGIN_AKANTU__

template<class Dumper>
class Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Dumpable(const std::string & filename) : dumper(filename) {};
  virtual ~Dumpable() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void addDumpField(const std::string & field_id) = 0;
  void removeDumpField(const std::string & field_id) {
    dumper.unRegisterField(field_id);
  };

  void setDirectory(const std::string & directory) {
    dumper.setDirectory(directory);
  }

  void setBaseName(const std::string & basename) {
    dumper.setBaseName(basename);
  }

  void dump() {
    dumper.dump();
  }


protected:
  void addDumpFieldToDumper(const std::string & field_id, DumperIOHelper::Field * field) {
    dumper.registerField(field_id, field);
  };

protected:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Dumper dumper;
};

#else

__BEGIN_AKANTU__
template<class Dumper>
class Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Dumpable(const std::string & filename) { };
  virtual ~Dumpable() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void addDumpField(const std::string & field_id) = 0;
  virtual void addDumpFieldVector(const std::string & field_id) {};
  virtual void addDumpFieldTensor(const std::string & field_id) {};
  void removeDumpField(const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  };

  void setDirectory(const std::string & directory) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void setBaseName(const std::string & basename) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

protected:
  void addDumpFieldToDumper(const std::string & field_id, Dumper::Field & field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  };

protected:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};
#endif //AKANTU_USE_IOHELPER

__END_AKANTU__

#endif /* __AKANTU_DUMPABLE_HH__ */
