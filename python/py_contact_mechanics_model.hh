#ifndef __AKANTU_PY_CONTACT_MECHANICS_MODEL_HH__
#define __AKANTU_PY_CONTACT_MECHANICS_MODEL_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_contact_mechanics_model(pybind11::module & mod);
} // namespace akantu

#endif //  __AKANTU_PY_CONTACT_MECHANICS_MODEL_HH__
