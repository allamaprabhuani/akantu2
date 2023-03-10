#define AKANTU_EIGEN_VERSION                                                   \
  (EIGEN_WORLD_VERSION * 10000 + EIGEN_MAJOR_VERSION * 1000 +                  \
   EIGEN_MAJOR_VERSION)

using size_type = Index;

template <typename PlainObjectType, int MapOptions, typename StrideType>
EIGEN_STRONG_INLINE
Matrix(const Map<PlainObjectType, MapOptions, StrideType> & other)
    : Base(other.derived().rows() * other.derived().cols(),
           other.derived().rows(), other.derived().cols()) {
  // AKANTU_DEBUG_WARNING("copy operator Map in matrix");
  Base::_check_template_params();
  Base::_resize_to_match(other);
  // FIXME/CHECK: isn't *this = other.derived() more efficient. it allows to
  //              go for pure _set() implementations, right?
  *this = other;
}

template <typename PlainObjectType, int MapOptions, typename StrideType>
EIGEN_STRONG_INLINE Matrix &
operator=(const Map<PlainObjectType, MapOptions, StrideType> & map) {
  // AKANTU_DEBUG_WARNING("operator= Map in matrix");
  return Base::_set(map);
}

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Matrix(std::initializer_list<Scalar> list)
    : Base(list) {}

#if AKANTU_EIGEN_VERSION <= 34000
template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
Matrix(std::initializer_list<std::initializer_list<Scalar>> list)
    : Base(list) {}
#endif
