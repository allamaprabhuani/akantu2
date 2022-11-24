#include "contact_detector_shared.hh"

namespace akantu {

AbstractContactDetector::AbstractContactDetector(akantu::Mesh & mesh)
    : mesh(mesh),
      spatial_dimension(mesh.getSpatialDimension()),
      positions(0, spatial_dimension) {

  this->positions.copy(mesh.getNodes());
}

} // namespace akantu
