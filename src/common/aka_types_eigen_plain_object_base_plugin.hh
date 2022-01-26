template <bool _is_vector = IsVectorAtCompileTime, typename _Derived = Derived,
          std::enable_if_t<_is_vector and
                           aka::is_eigen_matrix<_Derived>::value> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
PlainObjectBase(std::initializer_list<Scalar> list)
    : Base(list) {}

template <bool _is_vector = IsVectorAtCompileTime, typename _Derived = Derived,
          std::enable_if_t<not _is_vector and
                           aka::is_eigen_matrix<_Derived>::value> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
PlainObjectBase(std::initializer_list<std::initializer_list<Scalar>> list)
    : Base(list) {}
