#include "contact_mechanics_shared.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
AbstractContactDetector::AbstractContactDetector(akantu::Mesh & mesh, Array<Real> initial_positions)
    : mesh(mesh),
      spatial_dimension(mesh.getSpatialDimension()),
      positions(initial_positions) {
}

/* -------------------------------------------------------------------------- */
AbstractContactMechanicsModel::AbstractContactMechanicsModel(
    akantu::Mesh & mesh, const akantu::ModelType & type,
    std::shared_ptr<DOFManager> dof_manager, akantu::UInt dim,
    const akantu::ID id) : Model(mesh, type, dof_manager, dim, id) {
}

} // namespace akantu
