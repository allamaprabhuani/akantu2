#include "coupler_solid_contact_internodes.hh"

namespace akantu {

template <>
CouplerSolidContactInternodesTemplate<SolidMechanicsModel>::CouplerSolidContactInternodesTemplate(
    Mesh & mesh, UInt dim, const ID & id,
    std::shared_ptr<DOFManager> dof_manager)
    : AbstractCouplerSolidContactTemplate<SolidMechanicsModel,
                                          ContactMechanicsInternodesModel>(
          mesh, ModelType::_coupler_solid_contact, dim, id, dof_manager) {
  this->mesh.registerDumper<DumperParaview>("coupler_solid_contact", id, true);
  this->mesh.addDumpMeshToDumper("coupler_solid_contact", mesh,
                                 Model::spatial_dimension, _not_ghost,
                                 _ek_regular);
  this->registerDataAccessor(*this);

  solid = std::make_unique<SolidMechanicsModel>(mesh, Model::spatial_dimension,
                                                "solid_mechanics_model",
                                                this->dof_manager);
  contact = std::make_unique<ContactMechanicsInternodesModel>(
      mesh, Model::spatial_dimension, "contact_mechanics_internodes_model",
      this->dof_manager);
}

/* -------------------------------------------------------------------------- */
template <>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModel, ContactMechanicsInternodesModel>
    ::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  solid->initFull(_analysis_method = this->method);
  contact->initFull(_analysis_method = this->method);

  contact->setYoungsModulus(solid->getMaterial(0).get("E"));
}

} // namespace akantu
