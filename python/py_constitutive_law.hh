#include "py_aka_array.hh"
#include "py_akantu_pybind11_compatibility.hh"
/* -------------------------------------------------------------------------- */
#include <constitutive_law.hh>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PY_CONSTITUTIVE_LAW_HH_
#define AKANTU_PY_CONSTITUTIVE_LAW_HH_

namespace akantu {

void register_constitutive_law(pybind11::module & mod);

template <class ConstitutiveLawHandler>
class PyConstitutiveLaw : public ConstitutiveLaw<ConstitutiveLawHandler> {
public:
  using Parent = ConstitutiveLaw<ConstitutiveLawHandler>;
  /* Inherit the constructors */
  using Parent::Parent;

  ~PyConstitutiveLaw() override = default;

  void initConstitutiveLaw() override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, ConstitutiveLaw<ConstitutiveLawHandler>,
                      initConstitutiveLaw, );
  }

  template <typename T>
  void registerPyInternal(const std::string & name, UInt nb_component) {
    auto && internal = std::make_shared<InternalField<T>>(name, *this);
    AKANTU_DEBUG_INFO("alloc internal " << name << " " << internal);

    internal->initialize(nb_component);
    this->registerInternal(internal);
  }

  // methods need to be defined to be able to do the python interface
  UInt getNbData(const Array<Element> & elements,
                 const SynchronizationTag & tag) const override {
    return 0;
  }
  void packData(CommunicationBuffer & buffer,
                const Array<Element> & /*element*/,
                const SynchronizationTag & /*tag*/) const override {}
  void unpackData(CommunicationBuffer & /*buffer*/,
                  const Array<Element> & /*element*/,
                  const SynchronizationTag & /*tag*/) override{};
};

/* ------------------------------------------------------------------------ */
template <class ConstitutiveLawHandler>
void register_constitutive_law(py::module & mod, const std::string & name) {
  py::class_<ConstitutiveLaw<ConstitutiveLawHandler>,
             ConstitutiveLawInternalHandler, Parsable,
             PyConstitutiveLaw<ConstitutiveLawHandler>>(
      mod, name.c_str(), py::multiple_inheritance())
      .def(py::init<ConstitutiveLawHandler &, const ID &>())
      .def("registerInternalReal",
           [](ConstitutiveLaw<ConstitutiveLawHandler> & self,
              const std::string & name, UInt nb_component) {
             return dynamic_cast<PyConstitutiveLaw<ConstitutiveLawHandler> &>(
                        self)
                 .template registerPyInternal<Real>(name, nb_component);
           })
      .def("registerInternalUInt",
           [](ConstitutiveLaw<ConstitutiveLawHandler> & self,
              const std::string & name, UInt nb_component) {
             return dynamic_cast<PyConstitutiveLaw<ConstitutiveLawHandler> &>(
                        self)
                 .template registerPyInternal<UInt>(name, nb_component);
           })
      .def(
          "getInternalReal",
          [](ConstitutiveLaw<ConstitutiveLawHandler> & self, const ID & id)
              -> decltype(auto) { return self.template getInternal<Real>(id); },
          py::arg("id"), py::return_value_policy::reference)
      .def(
          "getInternalUInt",
          [](ConstitutiveLaw<ConstitutiveLawHandler> & self, const ID & id)
              -> decltype(auto) { return self.template getInternal<UInt>(id); },
          py::arg("id"), py::return_value_policy::reference)
      .def(
          "getElementFilter",
          [](ConstitutiveLaw<ConstitutiveLawHandler> & self) -> decltype(auto) {
            return self.getElementFilter();
          },
          py::return_value_policy::reference)
      .def("getSpatialDimension",
           &ConstitutiveLawInternalHandler::getSpatialDimension)
      /*
       * These functions override the `Parsable` interface.
       * This ensure that the `updateInternalParameters()` function is called.
       */
      .def(
          "setReal",
          [](ConstitutiveLaw<ConstitutiveLawHandler> & self, const ID & name,
             const Real value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setBool",
          [](ConstitutiveLaw<ConstitutiveLawHandler> & self, const ID & name,
             const bool value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setString",
          [](ConstitutiveLaw<ConstitutiveLawHandler> & self, const ID & name,
             const std::string & value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"))
      .def(
          "setInt",
          [](ConstitutiveLaw<ConstitutiveLawHandler> & self, const ID & name,
             const int value) -> void {
            self.setParam(name, value);
            return;
          },
          py::arg("name"), py::arg("value"));
}

} // namespace akantu

#endif // AKANTU_PY_CONSTITUTIVE_LAW_HH_
