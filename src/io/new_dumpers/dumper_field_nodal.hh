/* -------------------------------------------------------------------------- */
#include "dumper_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef DUMPER_FIELD_NODAL_H_
#define DUMPER_FIELD_NODAL_H_

namespace akantu {
namespace dumper {

  /* ------------------------------------------------------------------------ */
  class FieldNodalArrayBase : public FieldArrayBase {
  public:
    FieldNodalArrayBase(const SupportBase & support, std::type_index type)
        : FieldArrayBase(support, type, FieldType::_node_array) {}

    [[nodiscard]] Int size() const override {
      return aka::as_type<SupportElements>(this->getSupport()).getNbNodes();
    }
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Array_ = const Array<T> &>
  class FieldNodalArray : public FieldArray<T, Array_, FieldNodalArrayBase> {
  public:
    using FieldArray<T, Array_, FieldNodalArrayBase>::FieldArray;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class Function>
  class FieldFunctionNodalArray
      : public FieldFunctionArray<T, Function, FieldNodalArrayBase> {
  public:
    FieldFunctionNodalArray(
        const std::shared_ptr<FieldArrayTemplateBase<T, FieldNodalArrayBase>> &
            array_in,
        const SupportBase & support, Function && function) // NOLINT
        : FieldFunctionArray<T, Function, FieldNodalArrayBase>(
              array_in, support, std::forward<Function>(function)) {
      this->update();
    }

    FieldFunctionNodalArray(const Array<T> & array_in,
                            const SupportBase & support, Function && function)
        : FieldFunctionArray<T, Function, FieldNodalArrayBase>(
              array_in, support, std::forward<Function>(function)) {
      this->update();
    }
  };

  /* ------------------------------------------------------------------------ */
  template <class Array_, std::enable_if_t<dumper::details::is_array_v<
                              std::decay_t<Array_>>> * = nullptr>
  auto make_field(Array_ && array, const SupportBase & support) {
    using T = typename std::decay_t<Array_>::value_type;
    return std::make_shared<FieldNodalArray<T, Array_>>(
        std::forward<Array_>(array), support);
  }

  /* ------------------------------------------------------------------------ */
  template <class Function, class Array_,
            std::enable_if_t<
                dumper::details::is_array_v<std::decay_t<Array_>>> * = nullptr>
  auto make_field(Array_ && array, const SupportBase & support,
                  Function && function) {
    using T = typename std::decay_t<Array_>::value_type;
    return std::make_shared<FieldFunctionNodalArray<T, Function>>(
        std::forward<Array_>(array), support, std::forward<Function>(function));
  }

  /* ------------------------------------------------------------------------ */
  template <typename T, class Function>
  auto make_field(const std::shared_ptr<
                      FieldArrayTemplateBase<T, FieldNodalArrayBase>> & array,
                  const SupportBase & support, Function && function) {
    return std::make_shared<FieldFunctionNodalArray<T, Function>>(
        array, support, std::forward<Function>(function));
  }

} // namespace dumper
} // namespace akantu
#endif // DUMPER_FIELD_NODAL_H_
