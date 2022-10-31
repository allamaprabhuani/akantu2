/**
 * @file   test_tuples.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  mar déc 10 2019
 *
 * @brief A Documented file.
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * akantu-iterators is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * akantu-iterators is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with akantu-iterators. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_tuple_tools.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace aka;
using namespace aka::tuple;

TEST(NamedTuples, WithHash) {
  const auto a = std::vector<int>{1, 10, 3};
  const auto tuple = make_named_tuple(
      get<"nom"_h>() = std::string("Roger"), get<"age"_h>() = 47,
      get<"taille"_h>() = 1.92, get<"liste"_h>() = std::vector<int>({1, 2, 3}),
      get<"ref"_h>() = a);

  EXPECT_EQ(47, tuple.get(get<"age"_h>()));
  EXPECT_EQ("Roger", tuple.get(get<"nom"_h>()));
  EXPECT_EQ(1.92, tuple::get<"taille"_h>(tuple));
  EXPECT_EQ(a.data(), tuple::get<"ref"_h>(tuple).data());

  EXPECT_EQ(47, std::get<1>(tuple));
}

TEST(NamedTuples, StdGet) {
  const auto tuple = make_named_tuple(get<"a"_h>() = 1, get<"b"_h>() = 2);

  EXPECT_EQ(1, std::get<0>(tuple));
  EXPECT_EQ(2, std::get<1>(tuple));
}

TEST(NamedTuples, Transform) {
  const auto tuple = make_named_tuple(get<"a"_h>() = 1, get<"b"_h>() = 2);

  auto new_tuple = tuple::transform([](auto && a) { return 2 * a; }, tuple);

  EXPECT_EQ(2, tuple::get<"a"_h>(new_tuple));
  EXPECT_EQ(4, tuple::get<"b"_h>(new_tuple));
}

TEST(NamedTuples, Append) {
  const auto tuple = tuple::append(
      make_named_tuple(get<"a"_h>() = 1, get<"b"_h>() = 2), get<"c"_h>() = 3);

  EXPECT_EQ(1, tuple::get<"a"_h>(tuple));
  EXPECT_EQ(2, tuple::get<"b"_h>(tuple));
  EXPECT_EQ(3, tuple::get<"c"_h>(tuple));
}

TEST(NamedTuples, Replace) {
  const auto tuple = make_named_tuple(get<"a"_h>() = 1, get<"b"_h>() = 2);

  const auto new_tuple = tuple::replace<"a"_h>(tuple, 4);

  EXPECT_EQ(4, tuple::get<"a"_h>(new_tuple));
  EXPECT_EQ(2, tuple::get<"b"_h>(new_tuple));
}

TEST(NamedTuples, Has) {
  const auto tuple =
      tuple::append(make_named_tuple(get<"a"_h>() = 1, get<"b"_h>() = 2));

  auto has_a = tuple::has_t<"a"_h, decltype(tuple)>();
  auto has_b = tuple::has_t<"b"_h, decltype(tuple)>();
  auto has_c = tuple::has_t<"c"_h, decltype(tuple)>();

  static_assert(has_a, "Test of static assert on has_t");

  EXPECT_EQ(true, has_a);
  EXPECT_EQ(true, has_b);
  EXPECT_EQ(false, has_c);
}

#if defined(__INTEL_COMPILER)
// intel warnings here
#elif defined(__clang__)
// clang warnings here
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-string-literal-operator-template"
#elif (defined(__GNUC__) || defined(__GNUG__))
// gcc warnings here
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
TEST(NamedTuples, GNUExtension) {
  const auto a = std::vector<int>{1, 10, 3};
  const auto test = make_named_tuple(
      "nom"_n = std::string("Roger"), "age"_n = 47, "taille"_n = 1.92,
      "liste"_n = std::vector<int>({1, 2, 3}), "ref"_n = a);

  auto nom = test.get("nom"_n);
  EXPECT_EQ(47, test.get("age"_n));
  EXPECT_EQ("Roger", nom);
  EXPECT_EQ(1.92, test.get("taille"_n));
  EXPECT_EQ(a.data(), test.get("ref"_n).data());

  EXPECT_EQ(47, std::get<1>(test));
}

TEST(NamedTuples, MixGNUExtensionHash) {
  const auto tuple =
      make_named_tuple("nom"_n = std::string("Roger"), get<"age"_h>() = 47);

  EXPECT_EQ(47, tuple.get("age"_n));
}
#if defined(__clang__)
#pragma clang diagnostic pop
#elif (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic pop
#endif
