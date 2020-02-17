/* -------------------------------------------------------------------------- */
#include <fe_engine.hh>
#include <integration_point.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

__attribute__((visibility("default"))) void
register_fe_engine(py::module & mod) {

  py::class_<Element>(mod, "Element")
    .def(py::init<ElementType, UInt, GhostType>());

  py::class_<FEEngine>(mod, "FEEngine")
      .def("computeIntegrationPointsCoordinates",
           [](FEEngine & self, ElementTypeMapArray<Real> & coordinates,
              const ElementTypeMapArray<UInt> * filter_elements)
               -> decltype(auto) {
             return self.computeIntegrationPointsCoordinates(coordinates,
                                                             filter_elements);
           });

  py::class_<IntegrationPoint>(mod, "IntegrationPoint");
}
} // namespace akantu
