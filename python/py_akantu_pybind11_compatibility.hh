#ifndef PY_AKANTU_PYBIND11_COMPATIBILITY_HH_
#define PY_AKANTU_PYBIND11_COMPATIBILITY_HH_

#if not defined(PYBIND11_OVERRIDE)
#define PYBIND11_OVERRIDE PYBIND11_OVERLOAD
#define PYBIND11_OVERRIDE_NAME PYBIND11_OVERLOAD_NAME
#define PYBIND11_OVERRIDE_PURE PYBIND11_OVERLOAD_PURE
#define PYBIND11_OVERRIDE_PURE_NAME PYBIND11_OVERLOAD_PURE_NAME
#endif

#endif // PY_AKANTU_PYBIND11_COMPATIBILITY_HH_
