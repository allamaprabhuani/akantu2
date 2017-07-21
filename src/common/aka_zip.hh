/**
 * @file   aka_iterators.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Wed Jul 19 2017
 *
 * @brief iterator interfaces
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include <tuple>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_ITERATORS_HH__
#define __AKANTU_AKA_ITERATORS_HH__

namespace akantu {

namespace zip_ {
namespace details {
  struct dereference_iterator {
    template <class Iter> typename Iter::reference operator()(Iter & it) const {
      return *it;
    }
  };

  struct increment_iterator {
    template <class Iter> void operator()(Iter & it) const { ++it; }
  };

  struct begin_container {
    template <class Container>
    decltype(auto) operator()(Container && cont) const {
      return std::forward<Container>(cont).begin();
    }
  };

  struct end_container {
    template <class Container>
    decltype(auto) operator()(Container && cont) const {
      return std::forward<Container>(cont).end();
    }
  };
} // namespace details
} // namespace zip_

namespace tuple {
/* ------------------------------------------------------------------------ */
namespace details {
  template <size_t N> struct Foreach {
    template <class F, class Tuple>
    static inline decltype(auto) transform_forward(F && func, Tuple && tuple) {
      return std::tuple_cat(
          Foreach<N - 1>::transform_forward(std::forward<F>(func),
                                            std::forward<Tuple>(tuple)),
          std::forward_as_tuple(std::forward<F>(func)(
              std::get<N - 1>(std::forward<Tuple>(tuple)))));
    }

    template <class F, class Tuple>
    static inline decltype(auto) transform(F && func, Tuple && tuple) {
      return std::tuple_cat(
          Foreach<N - 1>::transform(std::forward<F>(func),
                                    std::forward<Tuple>(tuple)),
          std::make_tuple(std::forward<F>(func)(
              std::get<N - 1>(std::forward<Tuple>(tuple)))));
    }

    template <class F, class Tuple>
    static inline void foreach (F && func, Tuple && tuple) {
      Foreach<N - 1>::foreach (std::forward<F>(func),
                               std::forward<Tuple>(tuple));
      std::forward<F>(func)(std::get<N - 1>(std::forward<Tuple>(tuple)));
    }

    template <class Tuple> static inline bool equal(Tuple && a, Tuple && b) {
      if (not(std::get<N - 1>(std::forward<Tuple>(a)) ==
              std::get<N - 1>(std::forward<Tuple>(b))))
        return false;
      return Foreach<N - 1>::equal(std::forward<Tuple>(a),
                                   std::forward<Tuple>(b));
    }
  };

  /* ------------------------------------------------------------------------ */
  template <> struct Foreach<1> {
    template <class F, class Tuple>
    static inline decltype(auto) transform_forward(F && func, Tuple && tuple) {
      return std::forward_as_tuple(
          std::forward<F>(func)(std::get<0>(std::forward<Tuple>(tuple))));
    }

    template <class F, class Tuple>
    static inline decltype(auto) transform(F && func, Tuple && tuple) {
      return std::make_tuple(
          std::forward<F>(func)(std::get<0>(std::forward<Tuple>(tuple))));
    }

    template <class F, class Tuple>
    static inline void foreach (F && func, Tuple && tuple) {
      std::forward<F>(func)(std::get<0>(std::forward<Tuple>(tuple)));
    }

    template <class Tuple> static inline bool equal(Tuple && a, Tuple && b) {
      return std::get<0>(std::forward<Tuple>(a)) ==
             std::get<0>(std::forward<Tuple>(b));
    }
  };
} // namespace details
/* ------------------------------------------------------------------------ */
template <class Tuple> bool are_equal(Tuple && a, Tuple && b) {
  return details::Foreach<std::tuple_size<std::decay_t<Tuple>>::value>::equal(
      std::forward<Tuple>(a), std::forward<Tuple>(b));
}

template <class F, class Tuple> void foreach (F && func, Tuple && tuple) {
  details::Foreach<std::tuple_size<std::decay_t<Tuple>>::value>::foreach (
      std::forward<F>(func), std::forward<Tuple>(tuple));
}

template <class F, class Tuple>
decltype(auto) transform_forward(F && func, Tuple && tuple) {
  return details::Foreach<std::tuple_size<std::decay_t<Tuple>>::value>::
      transform_forward(std::forward<F>(func), std::forward<Tuple>(tuple));
}

template <class F, class Tuple>
decltype(auto) transform(F && func, Tuple && tuple) {
  return details::Foreach<std::tuple_size<std::decay_t<Tuple>>::value>::
      transform(std::forward<F>(func), std::forward<Tuple>(tuple));
}
} // namespace tuple

namespace zip_ {
/* -------------------------------------------------------------------------- */
template <class... Iterators> class iterator {
private:
  using tuple_t = std::tuple<Iterators...>;
public:
  explicit iterator(tuple_t iterators) : iterators(std::move(iterators)) {}

  decltype(auto) operator*() {
    return tuple::transform_forward(details::dereference_iterator(), iterators);
  }

  iterator & operator++() {
    tuple::foreach (details::increment_iterator(), iterators);
    return *this;
  }

  bool operator==(const iterator & other) const {
    return tuple::are_equal(iterators, other.iterators);
  }

  bool operator!=(const iterator & other) const {
    return not operator==(other);
  }

private:
  tuple_t iterators;
};

} // namespace zip_

template <class... Iterators>
decltype(auto) zip_iterator(std::tuple<Iterators...> && iterators_tuple) {
  auto zit = zip_::iterator<Iterators...>(std::move(iterators_tuple));
  return zit;
}

namespace zip_ {
template <class... Containers> class Zip {
  using containers_t = std::tuple<Containers &&...>;

public:
  explicit Zip(Containers &&... containers)
      : containers(std::forward_as_tuple(containers...)) {}

  decltype(auto) begin() const {
    return zip_iterator(tuple::transform(
        details::begin_container(), std::forward<containers_t>(containers)));
  }

  decltype(auto) end() const {
    return zip_iterator(tuple::transform(
        details::end_container(), std::forward<containers_t>(containers)));
  }

  decltype(auto) begin() {
    return zip_iterator(tuple::transform(
        details::begin_container(), std::forward<containers_t>(containers)));
  }

  decltype(auto) end() {
    return zip_iterator(tuple::transform(
        details::end_container(), std::forward<containers_t>(containers)));
  }

private:
  containers_t containers;
};
} // namespace zip_

template <class... Containers> decltype(auto) zip(Containers &&... containers) {
  return zip_::Zip<Containers...>(containers...);
}

} // namespace akantu

#endif /* __AKANTU_AKA_ITERATORS_HH__ */
