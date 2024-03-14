#include "aka_array.hh"
#include "element_type_map.hh"
#include "support.hh"
#include <type_traits>

#ifndef AKANTU_DUMPER_INTERNAL_TYPES_H_
#define AKANTU_DUMPER_INTERNAL_TYPES_H_

namespace akantu {
namespace dumper {

  namespace details {
    // primary template handles types that do not support pre-increment:
    template <class, class = void> constexpr bool has_getNbComponent_member{};

    // specialization recognizes types that do support pre-increment:
    template <class T>
    constexpr bool has_getNbComponent_member<
        T, std::void_t<decltype(std::declval<T &>().getNbComponent(1))>> = true;

    template <class T>
    constexpr bool has_getNbComponent_member<
        T, std::void_t<decltype(std::declval<T &>().getNbComponent(
               1, _not_defined, _casper))>> = true;

    /* ---------------------------------------------------------------------- */
    template <class, class = void> constexpr bool has_size_member{};

    template <class T>
    constexpr bool
        has_size_member<T, std::void_t<decltype(std::declval<T &>().size(1))>> =
            true;

    template <class T>
    constexpr bool
        has_size_member<T, std::void_t<decltype(std::declval<T &>().size(
                               1, _not_defined, _casper))>> = true;

    /* ---------------------------------------------------------------------- */
    template <class, class = void>
    constexpr bool has_set_nb_integration_points_member{};

    template <class T>
    constexpr bool has_set_nb_integration_points_member<
        T, std::void_t<decltype(std::declval<T &>().setNbIntegtrationPoints(
               ElementTypeMap<Int>{}))>> = true;

    /* ---------------------------------------------------------------------- */
    template <class T, class Function> struct function_with_type_return_scalar {
      using type = typename decltype(std::declval<Function &>().operator()(
          VectorProxy<T>(nullptr, 1), _not_defined, _casper))::Scalar;
    };

    template <class T, class Function>
    using function_with_type_return_scalar_t =
        typename function_with_type_return_scalar<T, Function>::type;

    /* ---------------------------------------------------------------------- */
    template <class T, class Function> struct function_return_scalar {
      using type = typename decltype(std::declval<Function &>().operator()(
          VectorProxy<T>(nullptr, 1)))::Scalar;
    };

    template <class T, class Function>
    using function_return_scalar_t =
        typename function_return_scalar<T, Function>::type;

    /* ---------------------------------------------------------------------- */
    template <class T> struct is_array : public std::false_type {};
    template <class T> struct is_array<Array<T>> : public std::true_type {};
    template <class T>
    struct is_array<const Array<T>> : public std::true_type {};

    template <class T> constexpr auto is_array_v = is_array<T>::value;

    template <class T>
    struct is_element_type_map_array : public std::false_type {};

    template <class T>
    struct is_element_type_map_array<ElementTypeMapArray<T>>
        : public std::true_type {};
    template <class T>
    struct is_element_type_map_array<const ElementTypeMapArray<T>>
        : public std::true_type {};
    template <class T>
    struct is_element_type_map_array<InternalField<T>> : public std::true_type {
    };

    template <class T>
    constexpr auto is_element_type_map_array_v =
        is_element_type_map_array<T>::value;

  } // namespace details

} // namespace dumper
} // namespace akantu

#endif // AKANTU_DUMPER_INTERNAL_TYPES_H_
