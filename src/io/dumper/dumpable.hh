/**
 * @file   dumpable.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 26 21:52:40 2012
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
#include "aka_common.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_paraview.hh"
#endif //AKANTU_USE_IOHELPER

#ifndef __AKANTU_DUMPABLE_HH__
#define __AKANTU_DUMPABLE_HH__

__BEGIN_AKANTU__

#ifdef AKANTU_USE_IOHELPER


template<class Dumper>
class Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Dumpable(const std::string & filename) : dumper(filename) { };
  virtual ~Dumpable() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void addDumpMesh(const Mesh & mesh, UInt spatial_dimension = _all_dimensions,
		   const GhostType & ghost_type = _not_ghost,
		   const ElementKind & element_kind = _ek_not_defined) {
    dumper.registerMesh(mesh,
			spatial_dimension,
			ghost_type,
			element_kind);
  }

  void addDumpFilteredMesh(const Mesh & mesh,
			   const ByElementTypeArray<UInt> & elements_filter,
			   const Array<UInt> & nodes_filter,		       
			   UInt spatial_dimension = _all_dimensions,
			   const GhostType & ghost_type = _not_ghost,
			   const ElementKind & element_kind = _ek_not_defined) {
    dumper.registerFilteredMesh(mesh,
				elements_filter,
				nodes_filter,
				spatial_dimension,
				ghost_type,
				element_kind);
  }

  virtual void addDumpField(const std::string & field_id) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  virtual void addDumpFieldExternal(__attribute__((unused)) const std::string & field_id,
				    __attribute__((unused)) DumperIOHelper::Field * field) {
    dumper.registerField(field_id, field);
  }

  template<typename T>
  void addDumpFieldExternal(const std::string & field_id, const Array<T> & field) {
    DumperIOHelper::Field * field_cont = new DumperIOHelper::NodalField<T>(field);
    dumper.registerField(field_id, field_cont);
  }

  template<typename T>
  void addDumpFieldExternal(const std::string & field_id,
			    const ByElementTypeArray<T> & field,
			    UInt spatial_dimension = _all_dimensions,
			    const GhostType & ghost_type = _not_ghost,
			    const ElementKind & element_kind = _ek_not_defined) {
    DumperIOHelper::Field * field_cont = new DumperIOHelper::ElementalField<T>(field,
									       spatial_dimension,
									       ghost_type,
									       element_kind);
    dumper.registerField(field_id, field_cont);
  }


  void removeDumpField(const std::string & field_id) {
    dumper.unRegisterField(field_id);
  }

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
  }

protected:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO_NOT_CONST(Dumper, dumper, Dumper &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Dumper dumper;
};

#else
/* -------------------------------------------------------------------------- */
class DumperIOHelper {

public:

  class Field;
};

class DumperParaview;
class Mesh;
class SubBoundary;
/* -------------------------------------------------------------------------- */


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

  void addDumpMesh(__attribute__((unused)) const Mesh & mesh,
		   __attribute__((unused)) UInt spatial_dimension = _all_dimensions,
		   __attribute__((unused)) const GhostType & ghost_type = _not_ghost,
		   __attribute__((unused)) const ElementKind & element_kind = _ek_not_defined) {
  }

  void addDumpFilteredMesh(__attribute__((unused)) const Mesh & mesh,
			   __attribute__((unused)) const SubBoundary & boundary,
			   __attribute__((unused)) UInt spatial_dimension = _all_dimensions,
			   __attribute__((unused)) const GhostType & ghost_type = _not_ghost,
			   __attribute__((unused)) const ElementKind & element_kind = _ek_not_defined) {
  }


  virtual void addDumpField(__attribute__((unused)) const std::string & field_id){
    AKANTU_DEBUG_TO_IMPLEMENT();
  };
  virtual void addDumpFieldExternal(__attribute__((unused)) const std::string & field_id,
				    __attribute__((unused)) DumperIOHelper::Field * field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  template<typename T>
  void addDumpFieldExternal(__attribute__((unused)) const std::string & field_id,
			    __attribute__((unused)) const Array<T> & field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  template<typename T>
  void addDumpFieldExternal(__attribute__((unused)) const std::string & field_id,
			    __attribute__((unused)) const ByElementTypeArray<T> & field,
			    __attribute__((unused)) UInt spatial_dimension = _all_dimensions,
			    __attribute__((unused)) const GhostType & ghost_type = _not_ghost,
			    __attribute__((unused)) const ElementKind & element_kind = _ek_not_defined) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  virtual void addDumpFieldArray(__attribute__((unused)) const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  };
  virtual void addDumpFieldTensor(__attribute__((unused)) const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  };
  void removeDumpField(__attribute__((unused)) const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void setDirectory(__attribute__((unused)) const std::string & directory) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void setBaseName(__attribute__((unused)) const std::string & basename) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void dump() {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

protected:

  void addDumpFieldToDumper(__attribute__((unused)) const std::string & field_id,
			    __attribute__((unused)) DumperIOHelper::Field * field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

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
