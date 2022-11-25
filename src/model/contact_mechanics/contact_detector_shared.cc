#include "contact_detector_shared.hh"

namespace akantu {

AbstractContactDetector::AbstractContactDetector(akantu::Mesh & mesh, Array<Real> initial_positions)
    : mesh(mesh),
      spatial_dimension(mesh.getSpatialDimension()),
      positions(initial_positions) {
}

} // namespace akantu
