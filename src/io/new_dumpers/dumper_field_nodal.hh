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
    FieldNodalArrayBase(SupportBase & support, std::type_index type)
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
  template <class T, class Function, class ReverseFunction>
  class FieldFunctionNodalArray
      : public FieldFunctionArray<T, Function, ReverseFunction,
                                  FieldNodalArrayBase> {
    using parent =
        FieldFunctionArray<T, Function, ReverseFunction, FieldNodalArrayBase>;

  public:
    FieldFunctionNodalArray(
        std::shared_ptr<FieldArrayTemplateBase<T, FieldNodalArrayBase>>
            array_in,
        SupportBase & support, Function && function,
        ReverseFunction && reverse_function) // NOLINT
        : parent(array_in, support, std::forward<Function>(function),
                 std::forward<ReverseFunction>(reverse_function)) {
      this->update();
    }

    FieldFunctionNodalArray(Array<T> & array_in, SupportBase & support,
                            Function && function,
                            ReverseFunction && reverse_function)
        : parent(array_in, support, std::forward<Function>(function),
                 std::forward<ReverseFunction>(reverse_function)) {
      this->update();
    }
  };

  /* ------------------------------------------------------------------------ */
  template <class Array_, std::enable_if_t<dumper::details::is_array_v<
                              std::decay_t<Array_>>> * = nullptr>
  auto make_field(Array_ && array, SupportBase & support) {
    using T = typename std::decay_t<Array_>::value_type;
    return std::make_shared<FieldNodalArray<T, Array_>>(
        std::forward<Array_>(array), support);
  }

  /* ------------------------------------------------------------------------ */
  template <
      typename T, class Function, class ReverseFunction = NoOpFunction,
      std::enable_if_t<std::is_invocable_v<Function, Vector<T>>> * = nullptr>
  auto make_field(
      std::shared_ptr<FieldArrayTemplateBase<T, FieldNodalArrayBase>> array,
      SupportBase & support, Function && function,
      ReverseFunction && reverse_function = NoOpFunction()) {
    return std::make_shared<
        FieldFunctionNodalArray<T, Function, ReverseFunction>>(
        array, support, std::forward<Function>(function),
        std::forward<ReverseFunction>(reverse_function));
  }

  /* ------------------------------------------------------------------------ */
  template <
      class Array_, class Function, class ReverseFunction = NoOpFunction,
      std::enable_if_t<
          dumper::details::is_array_v<std::decay_t<Array_>> and
          std::is_invocable_v<
              Function, Vector<typename std::decay_t<Array_>::value_type>>> * =
          nullptr>
  auto make_field(Array_ && array, SupportBase & support, Function && function,
                  ReverseFunction && reverse_function = NoOpFunction()) {
    using T = typename std::decay_t<Array_>::value_type;
    return std::make_shared<
        FieldFunctionNodalArray<T, Function, ReverseFunction>>(
        std::forward<Array_>(array), support, std::forward<Function>(function),
        std::forward<ReverseFunction>(reverse_function));
  }

} // namespace dumper
} // namespace akantu
#endif // DUMPER_FIELD_NODAL_H_
