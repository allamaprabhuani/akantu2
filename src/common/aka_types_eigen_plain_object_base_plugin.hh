
template <bool _is_vector = IsVectorAtCompileTime, typename _Derived = Derived,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC constexpr EIGEN_STRONG_INLINE
PlainObjectBase(std::initializer_list<Scalar> list) {
  static_assert(std::is_trivially_copyable<Scalar>{},
                "Cannot create a tensor on non trivial types");
  _check_template_params();
  this->template _init1<Index>(list.size());

  Index i = 0;
  for (auto val : list) {
    this->operator()(i++) = val;
  }
}

template <bool _is_vector = IsVectorAtCompileTime, typename _Derived = Derived,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC constexpr EIGEN_STRONG_INLINE
PlainObjectBase(std::initializer_list<std::initializer_list<Scalar>> list) {
  static_assert(std::is_trivially_copyable<Scalar>{},
                "Cannot create a tensor on non trivial types");
  Index m = list.size();
  Index n = 0;
  for (auto row : list) {
    n = std::max(n, Index(row.size()));
  }

  if (RowsAtCompileTime != -1 and RowsAtCompileTime != m) {
    throw std::range_error(
        "The size of the matrix does not correspond to the initializer_list");
  }

  if (ColsAtCompileTime != -1 and ColsAtCompileTime != n) {
    throw std::range_error(
        "The size of the matrix does not correspond to the initializer_list");
  }

  _check_template_params();
  this->template _init2<Index, Index>(m, n);
  this->fill(Scalar{});

  Index i = 0, j = 0;
  for (auto &row : list) {
    for (auto &val : row) {
      this->operator()(i, j++) = val;
    }
    ++i;
    j = 0;
  }
}
