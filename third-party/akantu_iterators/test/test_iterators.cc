/**
 * @file   test_zip_iterator.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 21 2017
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  test the zip container and iterator
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * akantu-iterators is free  software: you can redistribute it and/or  modify it
 * under the terms  of the  GNU Lesser  General Public  License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * akantu-iterators is  distributed in the  hope that it  will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public
 * License  for more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with akantu-iterators. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

using namespace aka;

/* -------------------------------------------------------------------------- */

// Non Trivial class that counts moves and copies
template <class T> class A {
public:
  A() = default;
  A(T a) : a(a){};
  A(const A & other)
      : a(other.a), copy_counter(other.copy_counter + 1),
        move_counter(other.move_counter) {}
  auto operator=(const A & other) -> A & {
    if (this != &other) {
      a = other.a;
      copy_counter = other.copy_counter + 1;
    }
    return *this;
  }

  A(A && other) noexcept
      : a(std::move(other.a)), copy_counter(std::move(other.copy_counter)),
        move_counter(std::move(other.move_counter) + 1) {}

  auto operator=(A && other) noexcept -> A & {
    if (this != &other) {
      a = std::move(other.a);
      copy_counter = std::move(other.copy_counter);
      move_counter = std::move(other.move_counter) + 1;
    }
    return *this;
  }

  auto operator*=(const T & b) -> A & {
    a *= b;
    return *this;
  }

  T a{};
  size_t copy_counter{0};
  size_t move_counter{0};
};

template <typename T> struct C {
  using size_type = std::size_t;
  struct iterator {
    using reference = A<T>;
    using difference_type = void;
    using iterator_category = std::input_iterator_tag;
    using value_type = A<T>;
    using pointer = A<T> *;

    iterator(T pos) : pos(std::move(pos)) {}

    auto operator*() -> A<T> { return A<int>(pos); }
    auto operator!=(const iterator & other) const -> bool {
      return pos != other.pos;
    }
    auto operator==(const iterator & other) const -> bool {
      return pos == other.pos;
    }
    auto operator++() -> iterator & {
      ++pos;
      return *this;
    }
    T pos;
  };

  C(T begin_, T end_) : begin_(std::move(begin_)), end_(std::move(end_)) {}

  auto begin() -> iterator { return iterator(begin_); }
  auto end() -> iterator { return iterator(end_); }

  T begin_, end_;
};

class TestZipFixutre : public ::testing::Test {
protected:
  void SetUp() override {
    a.reserve(size);
    b.reserve(size);

    for (size_t i = 0; i < size; ++i) {
      a.emplace_back(i);
      b.emplace_back(i + size);
    }
  }

  template <typename A, typename B>
  void check(A && a, B && b, size_t pos, size_t nb_copy, size_t nb_move) {
    EXPECT_EQ(pos, a.a);
    EXPECT_EQ(nb_copy, a.copy_counter);
    EXPECT_EQ(nb_move, a.move_counter);

    EXPECT_FLOAT_EQ(pos + this->size, b.a);
    EXPECT_EQ(nb_copy, b.copy_counter);
    EXPECT_EQ(nb_move, b.move_counter);
  }

public:
  std::size_t size{20};
  std::vector<A<int>> a{};
  std::vector<A<float>> b{};
};

TEST_F(TestZipFixutre, SimpleTest) {
  size_t i = 0;
  for (auto && pair : zip(this->a, this->b)) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 0);
    ++i;
  }
}

TEST_F(TestZipFixutre, ConstTest) {
  size_t i = 0;
  const auto & ca = this->a;
  const auto & cb = this->b;
  for (auto && pair : zip(ca, cb)) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 0);
    EXPECT_EQ(true,
              std::is_const<
                  std::remove_reference_t<decltype(std::get<0>(pair))>>::value);
    EXPECT_EQ(true,
              std::is_const<
                  std::remove_reference_t<decltype(std::get<1>(pair))>>::value);
    ++i;
  }
}

TEST_F(TestZipFixutre, MixteTest) {
  size_t i = 0;
  const auto & cb = this->b;
  for (auto && pair : zip(a, cb)) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 0);
    EXPECT_EQ(false,
              std::is_const<
                  std::remove_reference_t<decltype(std::get<0>(pair))>>::value);
    EXPECT_EQ(true,
              std::is_const<
                  std::remove_reference_t<decltype(std::get<1>(pair))>>::value);
    ++i;
  }
}

TEST_F(TestZipFixutre, MoveTest) {
  size_t i = 0;
  for (auto && pair :
       zip(C<int>(0, this->size), C<int>(this->size, 2 * this->size))) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 1);
    ++i;
  }
}

TEST_F(TestZipFixutre, Bidirectional) {
  auto _zip = zip(a, b);
  auto begin = _zip.begin();

  auto it = begin;
  ++it;
  EXPECT_EQ(begin, --it);

  it = begin;
  EXPECT_EQ(begin, it++);
  EXPECT_EQ(begin, --it);

  auto it2 = it = begin;
  ++it;
  ++it2;
  EXPECT_EQ(it2, it--);
  EXPECT_EQ(begin, it);
}

TEST_F(TestZipFixutre, RandomAccess) {
  auto _zip = zip(a, b);
  auto begin = _zip.begin();
  auto end = _zip.end();

  auto && val5 = begin[5];
  this->check(std::get<0>(val5), std::get<1>(val5), 5, 0, 0);

  auto && val13 = begin[13];
  this->check(std::get<0>(val13), std::get<1>(val13), 13, 0, 0);

  EXPECT_EQ(end - begin, a.size());
  auto it = ++begin;
  EXPECT_EQ(begin + 1, ++it);
  EXPECT_EQ(begin, it - 1);
}

TEST_F(TestZipFixutre, Cat) {
  size_t i = 0;
  for (auto && data : zip_cat(zip(a, b), zip(a, b))) {
    this->check(std::get<0>(data), std::get<1>(data), i, 0, 0);
    this->check(std::get<2>(data), std::get<3>(data), i, 0, 0);
    ++i;
  }
}

TEST_F(TestZipFixutre, AppendSimple) {
  size_t i = 0;

  auto zip_ = zip_append(zip(a, arange(size)), b);

  for (auto && data : zip_) {
    this->check(std::get<0>(data), std::get<2>(data), i, 0, 0);
    ++i;
  }
}

TEST_F(TestZipFixutre, AppendTemporaries) {
  size_t i = 0;

  auto zip_ = zip_append(zip_append(zip(a, arange(size)), b), arange(size));

  for (auto && data : zip_) {
    this->check(std::get<0>(data), std::get<2>(data), i, 0, 0);
    EXPECT_EQ(std::get<3>(data), i);
    ++i;
  }
}

TEST_F(TestZipFixutre, Replace) {
  size_t idx = 0;
  for (auto && [i, j] : zip_replace<1>(zip(a, a), b)) {
    this->check(i, j, idx, 0, 0);
    ++idx;
  }
}

TEST_F(TestZipFixutre, ReplaceReplace) {
  size_t i = 0;
  for (auto && data :
       zip_replace<1>(zip_replace<1>(zip(a, a), arange(size)), b)) {
    this->check(std::get<0>(data), std::get<1>(data), i, 0, 0);
    ++i;
  }
}

TEST_F(TestZipFixutre, Remove) {
  size_t i = 0;
  for (auto && data :
       zip_remove<1>(zip_remove<1>(zip(a, a, arange(size), b)))) {
    this->check(std::get<0>(data), std::get<1>(data), i, 0, 0);
    ++i;
  }
}

/* -------------------------------------------------------------------------- */
TEST_F(TestZipFixutre, SimpleNamedTest) {
  size_t i = 0;
  for (auto && pair : zip("a"_n = this->a, "b"_n = this->b)) {
    this->check(pair["a"_n], pair["b"_n], i, 0, 0);
    ++i;
  }
}

TEST_F(TestZipFixutre, NamedConstTest) {
  size_t i = 0;
  const auto & ca = this->a;
  const auto & cb = this->b;
  for (auto && pair : zip("a"_n = ca, "b"_n = cb)) {
    this->check(pair["a"_n], pair["b"_n], i, 0, 0);
    EXPECT_EQ(
        true,
        std::is_const<std::remove_reference_t<decltype(pair["a"_n])>>::value);
    EXPECT_EQ(
        true,
        std::is_const<std::remove_reference_t<decltype(pair["b"_n])>>::value);
    ++i;
  }
}

TEST_F(TestZipFixutre, NamedMixteTest) {
  size_t i = 0;
  const auto & cb = this->b;
  for (auto && pair : zip("a"_n = a, "b"_n = cb)) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 0);
    EXPECT_EQ(
        false,
        std::is_const<std::remove_reference_t<decltype(pair["a"_n])>>::value);
    EXPECT_EQ(
        true,
        std::is_const<std::remove_reference_t<decltype(pair["b"_n])>>::value);
    ++i;
  }
}

TEST_F(TestZipFixutre, MoveNamedTest) {
  size_t i = 0;
  for (auto && pair : zip("a"_n = C<int>(0, this->size),
                      "b"_n = C<int>(this->size, 2 * this->size))) {
    this->check(pair["a"_n], pair["b"_n], i, 0, 1);
    ++i;
  }
}

TEST_F(TestZipFixutre, BidirectionalNamed) {
  auto _zip = zip("a"_n = a, "b"_n = b);
  auto begin = _zip.begin();

  auto it = begin;
  ++it;
  EXPECT_EQ(begin, --it);

  it = begin;
  EXPECT_EQ(begin, it++);
  EXPECT_EQ(begin, --it);

  auto it2 = it = begin;
  ++it;
  ++it2;
  EXPECT_EQ(it2, it--);
  EXPECT_EQ(begin, it);
}

TEST_F(TestZipFixutre, RandomAccessNamed) {
  auto _zip = zip("a"_n = a, "b"_n = b);
  auto begin = _zip.begin();
  auto end = _zip.end();

  auto && val5 = begin[5];
  this->check(val5["a"_n], val5["b"_n], 5, 0, 0);

  auto && val13 = begin[13];
  this->check(val13["a"_n], val13["b"_n], 13, 0, 0);

  EXPECT_EQ(end - begin, a.size());
  auto it = ++begin;
  EXPECT_EQ(begin + 1, ++it);
  EXPECT_EQ(begin, it - 1);
}

TEST_F(TestZipFixutre, NamedAppend) {
  size_t i = 0;

  auto func = [&]() { return zip("a"_n = a, "temporary"_n = arange(size)); };

  auto && _zip = zip_append(zip_append(func(), "b"_n = b), "c"_n = arange(size),
                            "d"_n = b);
  for (auto && data : _zip) {
    this->check(data["a"_n], data["b"_n], i, 0, 0);
    EXPECT_EQ(data["c"_n], i);
    ++i;
  }
}

TEST_F(TestZipFixutre, NamedReplace) {
  size_t i = 0;
  for (auto && data :
       zip_replace(zip_replace(zip("a"_n = a, "b"_n = a, "c"_n = b),
                               "b"_n = arange(size)),
                   "b"_n = b)) {
    this->check(data["a"_n], data["b"_n], i, 0, 0);
    this->check(data["a"_n], data["c"_n], i, 0, 0);
    ++i;
  }
}

TEST_F(TestZipFixutre, NamedRemove) {
  size_t i = 0;
  for (auto && data :
       zip_remove(zip("a"_n = a, "a2"_n = arange(size), "b"_n = b), "a2"_n)) {
    EXPECT_EQ(data.has("a2"_n), false);
    this->check(data["a"_n], data["b"_n], i, 0, 0);
    ++i;
  }
}

TEST(TestNamedZipFixutre, Simple) {
  std::vector<int> a{0, 10, 20, 30, 40};
  std::vector<int> b{0, 1, 2, 3, 4};

  using namespace tuple;
  for (auto && data : zip("a"_n = a, "b"_n = b)) {
    auto & a = data["a"_n];
    auto & b = data["b"_n];
    b *= 10;
    EXPECT_EQ(b, a);
  }

  for (auto && data : zip("a"_n = a, "b"_n = b)) {
    auto & a = data["a"_n];
    auto & b = data["b"_n];
    EXPECT_EQ(b, a);
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestArangeIterator, Stop) {
  size_t ref_i = 0;
  for (auto i : arange(10)) {
    EXPECT_EQ(ref_i, i);
    ++ref_i;
  }
}

TEST(TestArangeIterator, StartStop) {
  size_t ref_i = 1;
  for (auto i : arange(1, 10)) {
    EXPECT_EQ(ref_i, i);
    ++ref_i;
  }
}

TEST(TestArangeIterator, StartStopStep) {
  size_t ref_i = 1;
  for (auto i : arange(1, 22, 2)) {
    EXPECT_EQ(ref_i, i);
    ref_i += 2;
  }
}

TEST(TestArangeIterator, StartStopStepZipped) {
  int ref_i1 = -1, ref_i2 = 1;
  for (auto && i : zip(arange(-1, -10, -1), arange(1, 18, 2))) {
    EXPECT_EQ(ref_i1, std::get<0>(i));
    EXPECT_EQ(ref_i2, std::get<1>(i));
    ref_i1 += -1;
    ref_i2 += 2;
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestEnumerateIterator, SimpleTest) {
  std::vector<int> a{0, 10, 20, 30, 40};
  std::vector<int> b{0, 2, 4, 6, 8};
  for (auto && data : enumerate(a, b)) {
    EXPECT_EQ(std::get<0>(data) * 10, std::get<1>(data));
    EXPECT_EQ(std::get<0>(data) * 2, std::get<2>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestTransformAdaptor, Keys) {
  std::map<std::string, int> map{
      {"1", 1}, {"2", 2}, {"3", 3}, {"3", 3}, {"4", 4}};

  char counter = '1';
  for (auto && key : make_keys_adaptor(map)) {
    EXPECT_EQ(counter, key[0]);
    ++counter;
  }
}

TEST(TestTransformAdaptor, Values) {
  std::map<std::string, int> map{
      {"1", 1}, {"2", 2}, {"3", 3}, {"3", 3}, {"4", 4}};

  int counter = 1;
  for (auto && value : make_values_adaptor(map)) {
    EXPECT_EQ(counter, value);
    ++counter;
  }
}

static int plus1(int value) { return value + 1; }

struct Plus {
  Plus(int a) : a(a) {}
  int operator()(int b) { return a + b; }

private:
  int a{0};
};

TEST(TestTransformAdaptor, Lambda) {
  auto && container = arange(10);

  for (auto && data :
       zip(container, make_transform_adaptor(container, [](auto && value) {
             return value + 1;
           }))) {
    EXPECT_EQ(std::get<0>(data) + 1, std::get<1>(data));
  }
}

TEST(TestTransformAdaptor, LambdaLambda) {
  std::map<std::string, int> map{
      {"1", 1}, {"2", 2}, {"3", 3}, {"3", 3}, {"4", 4}};

  int counter = 1;
  for (auto && data : make_transform_adaptor(
           make_values_adaptor(map), [](auto && value) { return value + 1; })) {
    EXPECT_EQ(counter + 1, data);
    ++counter;
  }

  auto && container = arange(10);

  for (auto && data :
       zip(container, make_transform_adaptor(container, [](auto && value) {
             return value + 1;
           }))) {
    EXPECT_EQ(std::get<0>(data) + 1, std::get<1>(data));
  }
}

TEST(TestTransformAdaptor, Function) {
  auto && container = arange(10);

  for (auto && data :
       zip(container, make_transform_adaptor(container, plus1))) {
    EXPECT_EQ(std::get<0>(data) + 1, std::get<1>(data));
  }
}

TEST(TestTransformAdaptor, Functor) {
  auto && container = arange(10);

  for (auto && data :
       zip(container, make_transform_adaptor(container, Plus(1)))) {
    EXPECT_EQ(std::get<0>(data) + 1, std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestFilteredIterator, Simple) {
  std::vector<int> values{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<int> filter_{1, 3, 4, 10, 8, 6};
  for (auto && data : zip(filter_, filter(filter_, values))) {
    EXPECT_EQ(std::get<0>(data), std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestFilteredIterator, Temporary) {
  std::vector<int> filter_{1, 3, 4, 10, 8, 6};
  for (auto && data :
       zip(filter_, filter(filter_, std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8,
                                                     9, 10}))) {
    EXPECT_EQ(std::get<0>(data), std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestFilteredIfIterator, Simple) {
  std::vector<int> values{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  for (auto && data : filter_if(values, [](auto && a) { return a % 2 == 0; })) {
    std::cout << data << std::endl;
    EXPECT_EQ(data % 2, 0);
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestFilteredIfIterator, Temporary) {
  for (auto && data :
       filter_if(arange(10), [](auto && a) { return a % 2 == 0; })) {
    EXPECT_EQ(data % 2, 0);
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestConcatenateIterator, SimpleTest) {
  for (auto && data : zip(arange(0, 13), concat(arange(0, 5), arange(5, 10),
                                                arange(10, 13)))) {
    EXPECT_EQ(std::get<0>(data), std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestRepeatNIterator, SimpleTest) {
  std::vector<int> ref{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4};
  for (auto && data : zip(repeat_n(arange(0, 5), 3), ref)) {
    EXPECT_EQ(std::get<0>(data), std::get<1>(data));
  }
}

TEST(TestRepeatNIterator, IncTest) {
  std::vector<int> ref{0, 1, 2, 3, 4};
  auto r = repeat_n(ref, 3);

  auto it = r.begin();
  EXPECT_EQ(*it, 0);
  EXPECT_EQ(it[3], 1);
  EXPECT_EQ(it[4], 1);
  EXPECT_EQ(it[5], 1);

  auto it1 = it + 1;
  EXPECT_EQ(it1[3], 1);
  EXPECT_EQ(it1[4], 1);
  EXPECT_EQ(it1[5], 2);

  auto end = r.end();
  EXPECT_EQ(*(it + 5), 1);
  EXPECT_EQ(*(end - 1), 4);
}
