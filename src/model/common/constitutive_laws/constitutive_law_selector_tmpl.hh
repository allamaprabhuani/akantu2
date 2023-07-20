#include "constitutive_law_selector.hh"

#ifndef AKANTU_CONSTITUTIVE_LAW_SELECTOR_TMPL_HH_
#define AKANTU_CONSTITUTIVE_LAW_SELECTOR_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, class Model_>
inline Idx ElementDataConstitutiveLawSelector<T, Model_>::operator()(
    const Element & element) {
  try {
    auto && data = this->elementData(element);
    if constexpr (std::is_same_v<T, std::string>) {
      return this->model.getConstitutiveLawIndex(data);
    } else {
      return data;
    }
  } catch (...) {
    return ConstitutiveLawSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, class Model_>
MeshDataConstitutiveLawSelector<T, Model_>::MeshDataConstitutiveLawSelector(
    const std::string & name, const Model_ & model, Idx first_index)
    : ElementDataConstitutiveLawSelector<T, Model_>(
          model.getMesh().template getData<T>(name), model, first_index) {}

} // namespace akantu

#endif // AKANTU_CONSTITUTIVE_LAW_SELECTOR_TMPL_HH_
