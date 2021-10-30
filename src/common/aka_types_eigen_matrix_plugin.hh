using size_type = Index;

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE
Matrix(std::initializer_list<std::initializer_list<Scalar>> list) {
  static_assert(std::is_trivially_copyable<Scalar>{},
                "Cannot create a tensor on non trivial types");

  StorageIndex n = 0;
  StorageIndex m = list.size();
  for (auto row : list) {
    n = std::max(n, row.size());
  }

  Base::_check_template_params();
  Base::template _init2<decltype(m), decltype(n)>(m, n);

  this->fill(Scalar{});

  decltype(m) i{0};
  decltype(n) j{0};
  for (auto & row : list) {
    for (auto & val : row) {
      (*this)[i, j++] = val;
    }
    ++i;
    j = 0;
  }
}

EIGEN_DEVICE_FUNC
EIGEN_STRONG_INLINE Matrix(std::initializer_list<Scalar> list) {
  static_assert(std::is_trivially_copyable<Scalar>{},
                "Cannot create a tensor on non trivial types");

  Base::_check_template_params();
  Base::template _init1<StorageIndex>(list.size());

  StorageIndex i = 0;
  for (auto val : list) {
    (*this)[i++] = val;
  }
}
