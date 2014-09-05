#include "dumpable.hh"

#ifdef AKANTU_USE_IOHELPER


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Dumpable::Dumpable() : default_dumper("") {
}

/* -------------------------------------------------------------------------- */
Dumpable::~Dumpable() {
}

/* -------------------------------------------------------------------------- */
void Dumpable::registerExternalDumper(DumperIOHelper & dumper,
                                             const std::string & dumper_name,
                                             const bool is_default) {
  this->dumpers[dumper_name] = &dumper;
  if (is_default)
    this->default_dumper = dumper_name;
}

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpMesh(const Mesh & mesh, UInt spatial_dimension,
                                  const GhostType & ghost_type,
                                  const ElementKind & element_kind) {

  this->addDumpMeshToDumper(this->default_dumper,
                            mesh,
                            spatial_dimension,
                            ghost_type,
                            element_kind);
}

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpMeshToDumper(const std::string & dumper_name,
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
void Dumpable::addDumpFilteredMesh(const Mesh & mesh,
                                          const ElementTypeMapArray<UInt> & elements_filter,
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
void Dumpable::addDumpFilteredMeshToDumper(const std::string & dumper_name,
                                                  const Mesh & mesh,
                                                  const ElementTypeMapArray<UInt> & elements_filter,
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
void Dumpable::addDumpField(const std::string & field_id) {
  this->addDumpFieldToDumper(this->default_dumper, field_id);
}

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpFieldToDumper(const std::string & dumper_name,
                                           const std::string & field_id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
};

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpFieldExternal(const std::string & field_id,
                                           dumper::Field * field) {
  this->addDumpFieldExternalToDumper(this->default_dumper, field_id, field);
}

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpFieldExternalToDumper(const std::string & dumper_name,
                                                   const std::string & field_id,
                                                   dumper::Field * field) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerField(field_id, field);
}

/* -------------------------------------------------------------------------- */

void Dumpable::removeDumpField(const std::string & field_id) {
  this->removeDumpFieldFromDumper(this->default_dumper, field_id);
}

/* -------------------------------------------------------------------------- */
void Dumpable::removeDumpFieldFromDumper(const std::string & dumper_name,
                                                const std::string & field_id) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.unRegisterField(field_id);
}

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpFieldVector(const std::string & field_id) {
  this->addDumpFieldVectorToDumper(this->default_dumper, field_id);
}

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpFieldVectorToDumper(const std::string & dumper_name,
                                                 const std::string & field_id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpFieldTensor(const std::string & field_id) {
  this->addDumpFieldTensorToDumper(this->default_dumper, field_id);
}

/* -------------------------------------------------------------------------- */
void Dumpable::addDumpFieldTensorToDumper(const std::string & dumper_name,
                                                 const std::string & field_id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void Dumpable::setDirectory(const std::string & directory) {
  this->setDirectoryToDumper(this->default_dumper, directory);
};

/* -------------------------------------------------------------------------- */
void Dumpable::setDirectoryToDumper(const std::string & dumper_name,
				    const std::string & directory) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.setDirectory(directory);
}

/* -------------------------------------------------------------------------- */
void Dumpable::setBaseName(const std::string & basename) {
  this->setBaseNameToDumper(this->default_dumper, basename);
}

/* -------------------------------------------------------------------------- */
void Dumpable::setBaseNameToDumper(const std::string & dumper_name,
                                          const std::string & basename) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.setBaseName(basename);
}

/* -------------------------------------------------------------------------- */
void Dumpable::setTimeStepToDumper(Real time_step) {
  this->setTimeStepToDumper(this->default_dumper, time_step);
}

/* -------------------------------------------------------------------------- */
void Dumpable::setTimeStepToDumper(const std::string & dumper_name,
                                          Real time_step) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.setTimeStep(time_step);
};

/* -------------------------------------------------------------------------- */
void Dumpable::dump(const std::string & dumper_name) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.dump();
}

/* -------------------------------------------------------------------------- */
void Dumpable::dump() {
  this->dump(this->default_dumper);
}

/* -------------------------------------------------------------------------- */
void Dumpable::dump(const std::string & dumper_name, UInt step) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.dump(step);
}

/* -------------------------------------------------------------------------- */
void Dumpable::dump(UInt step) {
  this->dump(this->default_dumper, step);
}

/* -------------------------------------------------------------------------- */
void Dumpable::dump(const std::string & dumper_name, Real time, UInt step) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.dump(time, step);
}

/* -------------------------------------------------------------------------- */
void Dumpable::dump(Real time, UInt step) {
  this->dump(this->default_dumper, time, step);
}


/* -------------------------------------------------------------------------- */
void Dumpable::internalAddDumpFieldToDumper(const std::string & dumper_name,
                                                   const std::string & field_id,
                                                   dumper::Field * field) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerField(field_id, field);
}


/* -------------------------------------------------------------------------- */
DumperIOHelper & Dumpable::getDumper() {
  return this->getDumper(this->default_dumper);
}


/* -------------------------------------------------------------------------- */
DumperIOHelper & Dumpable::getDumper(const std::string & dumper_name) {

  DumperMap::iterator it  = this->dumpers.find(dumper_name);
  DumperMap::iterator end = this->dumpers.end();

  if (it == end)
    AKANTU_EXCEPTION("Dumper " << dumper_name << "has not been registered, yet.");

  return *(it->second);
}

/* -------------------------------------------------------------------------- */

std::string Dumpable::getDefaultDumperName() const {
  return this->default_dumper;
}





__END_AKANTU__

#endif
