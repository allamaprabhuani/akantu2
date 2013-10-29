/**
 * @file   dumpable.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Kammer <david.kammer@epfl.ch>
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
#include "dumper_iohelper.hh"
#endif //AKANTU_USE_IOHELPER

#ifndef __AKANTU_DUMPABLE_HH__
#define __AKANTU_DUMPABLE_HH__

__BEGIN_AKANTU__

#ifdef AKANTU_USE_IOHELPER

class Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Dumpable() : default_dumper("") {};
  virtual ~Dumpable() {
    DumperMap::iterator it  = this->dumpers.begin();
    DumperMap::iterator end = this->dumpers.end();

    for (; it != end; ++it) {
      DumperSet::iterator fit  = this->external_dumpers.find(it->first);
      DumperSet::iterator fend = this->external_dumpers.end();

      if (fit == fend)
	delete it->second;
    }
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template<class T>
  void registerDumper(const std::string & dumper_name,
		      const std::string & file_name = "",
		      const bool is_default = false) {

    AKANTU_DEBUG_ASSERT(this->dumpers.find(dumper_name) == 
			this->dumpers.end(), 
			"Dumper " + dumper_name + "is already registered.");

    std::string name = file_name;
    if (name == "")
      name = dumper_name;

    this->dumpers[dumper_name] = new T(name);

    if (is_default)
      this->default_dumper = dumper_name;
  };

  void registerExternalDumper(DumperIOHelper * dumper,
			      const std::string & dumper_name,
			      const bool is_default = false) {
    this->dumpers[dumper_name] = dumper;
    if (is_default)
      this->default_dumper = dumper_name;
  };

  void addDumpMesh(const Mesh & mesh, UInt spatial_dimension = _all_dimensions,
		   const GhostType & ghost_type = _not_ghost,
		   const ElementKind & element_kind = _ek_not_defined) {

    this->addDumpMeshToDumper(this->default_dumper,
			      mesh,
			      spatial_dimension,
			      ghost_type,
			      element_kind);
  }

  void addDumpMeshToDumper(const std::string & dumper_name,
			   const Mesh & mesh, UInt spatial_dimension = _all_dimensions,
			   const GhostType & ghost_type = _not_ghost,
			   const ElementKind & element_kind = _ek_not_defined) {

    DumperIOHelper & dumper = this->getDumper(dumper_name);
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
    this->addDumpFilteredMeshToDumper(this->default_dumper,
				      mesh,
				      elements_filter,
				      nodes_filter,
				      spatial_dimension,
				      ghost_type,
				      element_kind);
  }

  void addDumpFilteredMeshToDumper(const std::string & dumper_name,
			   const Mesh & mesh,
			   const ByElementTypeArray<UInt> & elements_filter,
			   const Array<UInt> & nodes_filter,
			   UInt spatial_dimension = _all_dimensions,
			   const GhostType & ghost_type = _not_ghost,
			   const ElementKind & element_kind = _ek_not_defined) {

    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.registerFilteredMesh(mesh,
				elements_filter,
				nodes_filter,
				spatial_dimension,
				ghost_type,
				element_kind);
  };

  virtual void addDumpField(const std::string & field_id) {
    this->addDumpFieldToDumper(this->default_dumper, field_id);
  };

  virtual void addDumpFieldToDumper(const std::string & dumper_name, 
				    const std::string & field_id) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  virtual void addDumpFieldExternal(const std::string & field_id,
				    DumperIOHelper::Field * field) {
    this->addDumpFieldExternalToDumper(this->default_dumper, field_id, field);
  }

  virtual void addDumpFieldExternalToDumper(const std::string & dumper_name,
				    const std::string & field_id,
				    DumperIOHelper::Field * field) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.registerField(field_id, field);
  }

  template<typename T>
  void addDumpFieldExternal(const std::string & field_id, 
			    const Array<T> & field) {
    this->addDumpFieldExternalToDumper<T>(this->default_dumper,
					  field_id,
					  field);
  };
  
  template<typename T>
  void addDumpFieldExternalToDumper(const std::string & dumper_name,
				    const std::string & field_id, 
				    const Array<T> & field) {
    DumperIOHelper::Field * field_cont = new DumperIOHelper::NodalField<T>(field);
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.registerField(field_id, field_cont);
  }

  template<typename T>
  void addDumpFieldExternal(const std::string & field_id,
			    const ByElementTypeArray<T> & field,
			    UInt spatial_dimension = _all_dimensions,
			    const GhostType & ghost_type = _not_ghost,
			    const ElementKind & element_kind = _ek_not_defined) {
    this->addDumpFieldExternalToDumper(this->default_dumper,
				       field_id,
				       field,
				       spatial_dimension,
				       ghost_type,
				       element_kind);
  };

  template<typename T>
  void addDumpFieldExternalToDumper(const std::string & dumper_name,
				    const std::string & field_id,
				    const ByElementTypeArray<T> & field,
				    UInt spatial_dimension = _all_dimensions,
				    const GhostType & ghost_type = _not_ghost,
				    const ElementKind & element_kind = _ek_not_defined) {
    DumperIOHelper::Field * field_cont = new DumperIOHelper::ElementalField<T>(field,
									       spatial_dimension,
									       ghost_type,
									       element_kind);
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.registerField(field_id, field_cont);
  }


  void removeDumpField(const std::string & field_id) {
    this->removeDumpFieldFromDumper(this->default_dumper, field_id);
  }
  
  void removeDumpFieldFromDumper(const std::string & dumper_name,
			       const std::string & field_id) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.unRegisterField(field_id);
  }

  virtual void addDumpFieldVector(const std::string & field_id) {
    this->addDumpFieldVectorToDumper(this->default_dumper, field_id);
  }
  virtual void addDumpFieldVectorToDumper(const std::string & dumper_name,
					  const std::string & field_id) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  virtual void addDumpFieldTensor(const std::string & field_id) {
    this->addDumpFieldTensorToDumper(this->default_dumper, field_id);
  }
  virtual void addDumpFieldTensorToDumper(const std::string & dumper_name,
					  const std::string & field_id) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  void setDirectory(const std::string & directory) {
    this->setDirectoryToDumper(this->default_dumper, directory);
  };

  void setDirectoryToDumper(const std::string & dumper_name,
			    const std::string & directory) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.setDirectory(directory);
  }

  void setBaseName(const std::string & basename) {
    this->setBaseNameToDumper(this->default_dumper, basename);
  }

  void setBaseNameToDumper(const std::string & dumper_name,
			   const std::string & basename) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.setBaseName(basename);
  }

  void setTimeStepToDumper(Real time_step) {
    this->setTimeStepToDumper(this->default_dumper, time_step);
  }

  void setTimeStepToDumper(const std::string & dumper_name,
			   Real time_step) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.setTimeStep(time_step);
  };

  void dump(const std::string & dumper_name) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.dump();
  }

  void dump() {
    this->dump(this->default_dumper);
  }

  void dump(const std::string & dumper_name, UInt step) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.dump(step);
  }

  void dump(UInt step) {
    this->dump(this->default_dumper, step);
  }

  void dump(const std::string & dumper_name, Real time, UInt step) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.dump(time, step);
  }

  void dump(Real time, UInt step) {
    this->dump(this->default_dumper, time, step);
  }

protected:
  void internalAddDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id, 
				    DumperIOHelper::Field * field) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    dumper.registerField(field_id, field);
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  DumperIOHelper & getDumper() {
    return this->getDumper(this->default_dumper);
  }

  DumperIOHelper & getDumper(const std::string & dumper_name) {
    
    DumperMap::iterator it  = this->dumpers.find(dumper_name);
    DumperMap::iterator end = this->dumpers.end();
    
    if (it == end) 
      AKANTU_EXCEPTION("Dumper " << dumper_name << "has not been registered, yet.");

    return *(it->second);
  }

  template<class T>
  T & getDumper(const std::string & dumper_name) {
    DumperIOHelper & dumper = this->getDumper(dumper_name);
    
    try {
      T & templated_dumper = dynamic_cast<T & >(dumper); 
      return templated_dumper;
    }
    catch (...) {
      AKANTU_EXCEPTION("Dumper " << dumper_name << " is not of type: " 
		       << debug::demangle(typeid(T).name()));
    }
  }

  std::string getDefaultDumperName() const {
    return this->default_dumper;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::map<std::string, DumperIOHelper *> DumperMap;
  typedef std::set<std::string> DumperSet;

  DumperMap dumpers;
  std::string default_dumper;

  DumperSet external_dumpers;
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

class Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Dumpable() {};
  virtual ~Dumpable() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template<class T>
  void registerDumper(__attribute__((unused)) const std::string & dumper_name,
		      __attribute__((unused)) const std::string & file_name = "",
		      __attribute__((unused)) const bool is_default = false) {};

  void registerExternalDumper(__attribute__((unused)) DumperIOHelper * dumper,
			      __attribute__((unused)) const std::string & dumper_name,
			      __attribute__((unused)) const bool is_default = false) {
  };

  void addDumpMesh(__attribute__((unused)) const Mesh & mesh,
		   __attribute__((unused)) UInt spatial_dimension = _all_dimensions,
		   __attribute__((unused)) const GhostType & ghost_type = _not_ghost,
		   __attribute__((unused)) const ElementKind & element_kind = _ek_not_defined) {
  }

  void addDumpMeshToDumper(__attribute__((unused)) const std::string & dumper_name,
			   __attribute__((unused)) const Mesh & mesh,
			   __attribute__((unused)) UInt spatial_dimension = _all_dimensions,
			   __attribute__((unused)) const GhostType & ghost_type = _not_ghost,
			   __attribute__((unused)) const ElementKind & element_kind = _ek_not_defined) {
  }

  void addDumpFilteredMesh(__attribute__((unused)) const Mesh & mesh,
			   __attribute__((unused)) const ByElementTypeArray<UInt> & elements_filter,
			   __attribute__((unused)) const Array<UInt> & nodes_filter,
			   __attribute__((unused)) UInt spatial_dimension = _all_dimensions,
			   __attribute__((unused)) const GhostType & ghost_type = _not_ghost,
			   __attribute__((unused)) const ElementKind & element_kind = _ek_not_defined) {
  }

  void addDumpFilteredMeshToDumper(__attribute__((unused)) const std::string & dumper_name,
				   __attribute__((unused)) const Mesh & mesh,
				   __attribute__((unused)) const ByElementTypeArray<UInt> & elements_filter,
				   __attribute__((unused)) const Array<UInt> & nodes_filter,
				   __attribute__((unused)) UInt spatial_dimension = _all_dimensions,
				   __attribute__((unused)) const GhostType & ghost_type = _not_ghost,
				   __attribute__((unused)) const ElementKind & element_kind = _ek_not_defined) {
  }

  virtual void addDumpField(__attribute__((unused)) const std::string & field_id){
    AKANTU_DEBUG_TO_IMPLEMENT();
  };
  virtual void addDumpFieldToDumper(__attribute__((unused)) const std::string & dumper_name,
				    __attribute__((unused)) const std::string & field_id){
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  virtual void addDumpFieldExternal(__attribute__((unused)) const std::string & field_id,
				    __attribute__((unused)) DumperIOHelper::Field * field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  virtual void addDumpFieldExternalToDumper(__attribute__((unused)) const std::string & dumper_name,
				    __attribute__((unused)) const std::string & field_id,
				    __attribute__((unused)) DumperIOHelper::Field * field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  template<typename T>
  void addDumpFieldExternal(__attribute__((unused)) const std::string & field_id,
			    __attribute__((unused)) const Array<T> & field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  template<typename T>
  void addDumpFieldExternalToDumper(__attribute__((unused)) const std::string & dumper_name,
				    __attribute__((unused)) const std::string & field_id,
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
  template<typename T>
  void addDumpFieldExternalToDumper(__attribute__((unused)) const std::string & dumper_name,
				    __attribute__((unused)) const std::string & field_id,
				    __attribute__((unused)) const ByElementTypeArray<T> & field,
				    __attribute__((unused)) UInt spatial_dimension = _all_dimensions,
				    __attribute__((unused)) const GhostType & ghost_type = _not_ghost,
				    __attribute__((unused)) const ElementKind & element_kind = _ek_not_defined) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void removeDumpField(__attribute__((unused)) const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void removeDumpFieldFromDumper(__attribute__((unused)) const std::string & dumper_name,
			       __attribute__((unused)) const std::string & field_id) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void setDirectory(__attribute__((unused)) const std::string & directory) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void setDirectoryToDumper(__attribute__((unused)) const std::string & dumper_name,
			    __attribute__((unused)) const std::string & directory) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void setBaseName(__attribute__((unused)) const std::string & basename) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void setBaseNameToDumper(__attribute__((unused)) const std::string & dumper_name,
			   __attribute__((unused)) const std::string & basename) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void dump() {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(__attribute__((unused)) const std::string & dumper_name) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void dump(__attribute__((unused)) UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(__attribute__((unused)) const std::string & dumper_name,
	    __attribute__((unused)) UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  void dump(__attribute__((unused)) Real current_time,
	    __attribute__((unused)) UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }
  void dump(__attribute__((unused)) const std::string & dumper_name,
	    __attribute__((unused)) Real current_time,
	    __attribute__((unused)) UInt step) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

protected:

  void internalAddDumpFieldToDumper(__attribute__((unused)) const std::string & dumper_name,
				    __attribute__((unused)) const std::string & field_id,
				    __attribute__((unused)) DumperIOHelper::Field * field) {
    AKANTU_DEBUG_WARNING("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

protected:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  DumperIOHelper & getDumper() {
    AKANTU_DEBUG_ERROR("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  }

  DumperIOHelper & getDumper(__attribute__((unused)) const std::string & dumper_name){
    AKANTU_DEBUG_ERROR("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  };

  template<class T>
  T & getDumper(__attribute__((unused)) const std::string & dumper_name) {
    AKANTU_DEBUG_ERROR("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  };

  std::string getDefaultDumperName() {
    AKANTU_DEBUG_ERROR("No dumper activated at compilation, turn on AKANTU_USE_IOHELPER in cmake.");
  };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};
#endif //AKANTU_USE_IOHELPER

__END_AKANTU__

#endif /* __AKANTU_DUMPABLE_HH__ */
