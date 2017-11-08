/**
 * @file   dumper_iohelper.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 26 2012
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  Define the akantu dumper interface for IOhelper dumpers
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "aka_types.hh"
#include "aka_array.hh"
#include "element_type_map.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DUMPER_IOHELPER_HH__
#define __AKANTU_DUMPER_IOHELPER_HH__
/* -------------------------------------------------------------------------- */

namespace iohelper {
class Dumper;
}

namespace akantu {

UInt getIOHelperType(ElementType type);

namespace dumper {
class Field;
class VariableBase;
}

class Mesh;

class DumperIOHelper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DumperIOHelper();
  virtual ~DumperIOHelper();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register a given Mesh for the current dumper
  virtual void registerMesh(const Mesh & mesh,
                            UInt spatial_dimension = _all_dimensions,
                            const GhostType & ghost_type = _not_ghost,
                            const ElementKind & element_kind = _ek_not_defined);

  /// register a filtered Mesh (provided filter lists) for the current dumper
  virtual void
  registerFilteredMesh(const Mesh & mesh,
                       const ElementTypeMapArray<UInt> & elements_filter,
                       const Array<UInt> & nodes_filter,
                       UInt spatial_dimension = _all_dimensions,
                       const GhostType & ghost_type = _not_ghost,
                       const ElementKind & element_kind = _ek_not_defined);

  /// register a Field object identified by name and provided by pointer
  void registerField(const std::string & field_id, dumper::Field * field);
  /// remove the Field identified by name from managed fields
  void unRegisterField(const std::string & field_id);
  /// register a VariableBase object identified by name and provided by pointer
  void registerVariable(const std::string & variable_id,
                        dumper::VariableBase * variable);
  /// remove a VariableBase identified by name from managed fields
  void unRegisterVariable(const std::string & variable_id);

  /// request dump: this calls IOHelper dump routine
  virtual void dump();
  /// request dump: this first set the current step and then calls IOHelper dump
  /// routine
  virtual void dump(UInt step);
  /// request dump: this first set the current step and current time and then
  /// calls IOHelper dump routine
  virtual void dump(Real current_time, UInt step);
  /// set the parallel context for IOHeper
  virtual void setParallelContext(bool is_parallel);
  /// set the directory where to generate the dumped files
  virtual void setDirectory(const std::string & directory);
  /// set the base name (needed by most IOHelper dumpers)
  virtual void setBaseName(const std::string & basename);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// direct access to the iohelper::Dumper object
  AKANTU_GET_MACRO(Dumper, *dumper, iohelper::Dumper &)

  /// set the timestep of the iohelper::Dumper
  void setTimeStep(Real time_step);

public:
  /* ------------------------------------------------------------------------ */
  /* Variable wrapper */
  template <typename T, bool is_scal = std::is_arithmetic<T>::value> class Variable;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// internal iohelper::Dumper
  iohelper::Dumper * dumper;

  typedef std::map<std::string, dumper::Field *> Fields;
  typedef std::map<std::string, dumper::VariableBase *> Variables;

  /// list of registered fields to dump
  Fields fields;
  Variables variables;

  /// dump counter
  UInt count{0};

  /// directory name
  std::string directory;

  /// filename prefix
  std::string filename;

  /// is time tracking activated in the dumper
  bool time_activated{false};
};

} // akantu

#endif /* __AKANTU_DUMPER_IOHELPER_HH__ */
