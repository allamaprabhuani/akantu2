using size_type = Index;

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
Matrix(std::initializer_list<std::initializer_list<Scalar>> list)
    : Base(list) {}

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Matrix(std::initializer_list<Scalar> list)
    : Base(list) {}
