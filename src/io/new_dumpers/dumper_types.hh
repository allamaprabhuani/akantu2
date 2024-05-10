/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_enum_macros.hh"
/* -------------------------------------------------------------------------- */
#include <cstdint>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef DUMPER_TYPES_HH
#define DUMPER_TYPES_HH

namespace akantu {
namespace dumper {
  class FieldBase;
  class FieldNodalArrayBase;
  class FieldElementalArrayBase;
  template <class T, class Base> class FieldArrayTemplateBase;
  template <class T>
  using FieldNodalArrayTemplateBase =
      FieldArrayTemplateBase<T, FieldNodalArrayBase>;
  template <class T>
  using FieldElementalArrayTemplateBase =
      FieldArrayTemplateBase<T, FieldElementalArrayBase>;
  template <class T> class FieldElementMapArrayTemplateBase;

  class VariableBase;

  // clang-format off
  enum class FieldUsageType : uint8_t {
    _internal      = 0b10000000,
    _restart       = 0b00000001,
    _visualisation = 0b00000010
  };
  // clang-format on

  /// Bit-wise operator between access modes
  inline FieldUsageType operator|(const FieldUsageType & a,
                                  const FieldUsageType & b) {
    using ut = std::underlying_type_t<FieldUsageType>;
    auto tmp = FieldUsageType(ut(a) | ut(b));
    return tmp;
  }

  inline bool operator&(const FieldUsageType & a, const FieldUsageType & b) {
    using ut = std::underlying_type_t<FieldUsageType>;
    auto tmp = FieldUsageType(ut(a) & ut(b));
    return ut(tmp) != 0;
  }

#define AKANTU_SUPPORT_TYPES (mesh)(element_group)

  AKANTU_CLASS_ENUM_DECLARE(SupportType, AKANTU_SUPPORT_TYPES)

  template <SupportType t>
  using support_type_t = std::integral_constant<SupportType, t>;

  using support_type_mesh_t = support_type_t<SupportType::_mesh>;
  using support_type_element_group_t =
      support_type_t<SupportType::_element_group>;

  using AllSupportTypes =
      std::tuple<support_type_mesh_t, support_type_element_group_t>;

#define AKANTU_FIELD_TYPES                                                     \
  (not_defined)(array)(                                                        \
      node_array)(element_array)(element_map_array)(internal_field)

  AKANTU_CLASS_ENUM_DECLARE(FieldType, AKANTU_FIELD_TYPES)

  template <FieldType t>
  using field_type_t = std::integral_constant<FieldType, t>;

// creating a type instead of a using helps to debug
#define AKANTU_DECLARE_FIELD_TYPES(r, data, ty)                                \
  using BOOST_PP_CAT(BOOST_PP_CAT(field_type_, ty), _t) =                      \
      field_type_t<BOOST_PP_CAT(FieldType::_, ty)>;
  BOOST_PP_SEQ_FOR_EACH(AKANTU_DECLARE_FIELD_TYPES, _, AKANTU_FIELD_TYPES)

#define OP_CAT(s, data, ty) field_type_t<BOOST_PP_CAT(FieldType::_, ty)>
  using AllFieldTypes = std::tuple<BOOST_PP_SEQ_ENUM(
      BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_FIELD_TYPES))>;
#undef OP_CAT
} // namespace dumper

AKANTU_CLASS_ENUM_INPUT_STREAM(dumper::FieldType, AKANTU_FIELD_TYPES)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(dumper::FieldType, AKANTU_FIELD_TYPES)
AKANTU_CLASS_ENUM_INPUT_STREAM(dumper::SupportType, AKANTU_SUPPORT_TYPES)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(dumper::SupportType, AKANTU_SUPPORT_TYPES)

} // namespace akantu

#endif // DUMPER_TYPES_HH
