#ifndef __AKANTU_PY_MODEL_COUPLERS_HH__
#define __AKANTU_PY_MODEL_COUPLERS_HH__

namespace pybind11 {
struct module;
} // namespace pybind11

namespace akantu {
void register_model_couplers(pybind11::module & mod);
} // namespace akantu

#endif //  __AKANTU_PY_MODEL_COUPLERS_HH__
