/**
 * @file   dumpable_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Sat Jul 13 11:35:22 2013
 *
 * @brief  Implementation of the Dumpable class
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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline Dumpable::Dumpable() : default_dumper("") {
}

/* -------------------------------------------------------------------------- */
inline Dumpable::~Dumpable() {
  DumperMap::iterator it  = this->dumpers.begin();
  DumperMap::iterator end = this->dumpers.end();

  for (; it != end; ++it) {
    DumperSet::iterator fit  = this->external_dumpers.find(it->first);
    DumperSet::iterator fend = this->external_dumpers.end();

    if (fit == fend)
      delete it->second;
  }
}

/* -------------------------------------------------------------------------- */
template<class T>
inline void Dumpable::registerDumper(const std::string & dumper_name,
                                     const std::string & file_name,
                                     const bool is_default) {

  AKANTU_DEBUG_ASSERT(this->dumpers.find(dumper_name) ==
                      this->dumpers.end(),
                      "Dumper " + dumper_name + "is already registered.");

  std::string name = file_name;
  if (name == "")
    name = dumper_name;

  this->dumpers[dumper_name] = new T(name);

  if (is_default)
    this->default_dumper = dumper_name;
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::registerExternalDumper(DumperIOHelper & dumper,
                                             const std::string & dumper_name,
                                             const bool is_default) {
  this->dumpers[dumper_name] = &dumper;
  if (is_default)
    this->default_dumper = dumper_name;
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpMesh(const Mesh & mesh, UInt spatial_dimension,
                                  const GhostType & ghost_type,
                                  const ElementKind & element_kind) {

  this->addDumpMeshToDumper(this->default_dumper,
                            mesh,
                            spatial_dimension,
                            ghost_type,
                            element_kind);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpMeshToDumper(const std::string & dumper_name,
                                          const Mesh & mesh, UInt spatial_dimension,
                                          const GhostType & ghost_type,
                                          const ElementKind & element_kind) {

  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerMesh(mesh,
                      spatial_dimension,
                      ghost_type,
                      element_kind);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFilteredMesh(const Mesh & mesh,
                                          const ByElementTypeArray<UInt> & elements_filter,
                                          const Array<UInt> & nodes_filter,
                                          UInt spatial_dimension,
                                          const GhostType & ghost_type,
                                          const ElementKind & element_kind) {
  this->addDumpFilteredMeshToDumper(this->default_dumper,
                                    mesh,
                                    elements_filter,
                                    nodes_filter,
                                    spatial_dimension,
                                    ghost_type,
                                    element_kind);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFilteredMeshToDumper(const std::string & dumper_name,
                                                  const Mesh & mesh,
                                                  const ByElementTypeArray<UInt> & elements_filter,
                                                  const Array<UInt> & nodes_filter,
                                                  UInt spatial_dimension,
                                                  const GhostType & ghost_type,
                                                  const ElementKind & element_kind) {

  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerFilteredMesh(mesh,
                              elements_filter,
                              nodes_filter,
                              spatial_dimension,
                              ghost_type,
                              element_kind);
}


/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpField(const std::string & field_id) {
  this->addDumpFieldToDumper(this->default_dumper, field_id);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFieldToDumper(const std::string & dumper_name,
                                           const std::string & field_id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
};

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFieldExternal(const std::string & field_id,
                                           DumperIOHelper::Field * field) {
  this->addDumpFieldExternalToDumper(this->default_dumper, field_id, field);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFieldExternalToDumper(const std::string & dumper_name,
                                                   const std::string & field_id,
                                                   DumperIOHelper::Field * field) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerField(field_id, field);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Dumpable::addDumpFieldExternal(const std::string & field_id,
                                           const Array<T> & field) {
  this->addDumpFieldExternalToDumper<T>(this->default_dumper,
                                        field_id,
                                        field);
};

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Dumpable::addDumpFieldExternalToDumper(const std::string & dumper_name,
                                                   const std::string & field_id,
                                                   const Array<T> & field) {
  DumperIOHelper::Field * field_cont = new DumperIOHelper::NodalField<T>(field);
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerField(field_id, field_cont);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Dumpable::addDumpFieldExternal(const std::string & field_id,
                                           const ByElementTypeArray<T> & field,
                                           UInt spatial_dimension,
                                           const GhostType & ghost_type,
                                           const ElementKind & element_kind) {
  this->addDumpFieldExternalToDumper(this->default_dumper,
                                     field_id,
                                     field,
                                     spatial_dimension,
                                     ghost_type,
                                     element_kind);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Dumpable::addDumpFieldExternalToDumper(const std::string & dumper_name,
                                                   const std::string & field_id,
                                                   const ByElementTypeArray<T> & field,
                                                   UInt spatial_dimension,
                                                   const GhostType & ghost_type,
                                                   const ElementKind & element_kind) {
  DumperIOHelper::Field * field_cont = new DumperIOHelper::ElementalField<T>(field,
                                                                             spatial_dimension,
                                                                             ghost_type,
                                                                             element_kind);
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerField(field_id, field_cont);
}


/* -------------------------------------------------------------------------- */
inline void Dumpable::removeDumpField(const std::string & field_id) {
  this->removeDumpFieldFromDumper(this->default_dumper, field_id);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::removeDumpFieldFromDumper(const std::string & dumper_name,
                                                const std::string & field_id) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.unRegisterField(field_id);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFieldVector(const std::string & field_id) {
  this->addDumpFieldVectorToDumper(this->default_dumper, field_id);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFieldVectorToDumper(const std::string & dumper_name,
                                                 const std::string & field_id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFieldTensor(const std::string & field_id) {
  this->addDumpFieldTensorToDumper(this->default_dumper, field_id);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::addDumpFieldTensorToDumper(const std::string & dumper_name,
                                                 const std::string & field_id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::setDirectory(const std::string & directory) {
  this->setDirectoryToDumper(this->default_dumper, directory);
};

/* -------------------------------------------------------------------------- */
inline void Dumpable::setDirectoryToDumper(const std::string & dumper_name,
                                           const std::string & directory) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.setDirectory(directory);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::setBaseName(const std::string & basename) {
  this->setBaseNameToDumper(this->default_dumper, basename);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::setBaseNameToDumper(const std::string & dumper_name,
                                          const std::string & basename) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.setBaseName(basename);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::setTimeStepToDumper(Real time_step) {
  this->setTimeStepToDumper(this->default_dumper, time_step);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::setTimeStepToDumper(const std::string & dumper_name,
                                          Real time_step) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.setTimeStep(time_step);
};

/* -------------------------------------------------------------------------- */
inline void Dumpable::dump(const std::string & dumper_name) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.dump();
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::dump() {
  this->dump(this->default_dumper);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::dump(const std::string & dumper_name, UInt step) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.dump(step);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::dump(UInt step) {
  this->dump(this->default_dumper, step);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::dump(const std::string & dumper_name, Real time, UInt step) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.dump(time, step);
}

/* -------------------------------------------------------------------------- */
inline void Dumpable::dump(Real time, UInt step) {
  this->dump(this->default_dumper, time, step);
}


/* -------------------------------------------------------------------------- */
inline void Dumpable::internalAddDumpFieldToDumper(const std::string & dumper_name,
                                                   const std::string & field_id,
                                                   DumperIOHelper::Field * field) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerField(field_id, field);
}


/* -------------------------------------------------------------------------- */
inline DumperIOHelper & Dumpable::getDumper() {
  return this->getDumper(this->default_dumper);
}


/* -------------------------------------------------------------------------- */
inline DumperIOHelper & Dumpable::getDumper(const std::string & dumper_name) {

  DumperMap::iterator it  = this->dumpers.find(dumper_name);
  DumperMap::iterator end = this->dumpers.end();

  if (it == end)
    AKANTU_EXCEPTION("Dumper " << dumper_name << "has not been registered, yet.");

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template<class T>
inline T & Dumpable::getDumper(const std::string & dumper_name) {
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

/* -------------------------------------------------------------------------- */
inline std::string Dumpable::getDefaultDumperName() const {
  return this->default_dumper;
}


__END_AKANTU__
