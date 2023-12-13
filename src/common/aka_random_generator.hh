/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <random>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_RANDOM_GENERATOR_HH_
#define AKANTU_AKA_RANDOM_GENERATOR_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* List of available distributions                                            */
/* -------------------------------------------------------------------------- */
#define AKANTU_RANDOM_DISTRIBUTION_TYPES                                       \
  (uniform)(exponential)(                                                      \
      gamma)(weibull)(extreme_value)(normal)(lognormal)(chi_squared)(cauchy)(fisher_f)(student_t)(not_defined)

AKANTU_CLASS_ENUM_DECLARE(RandomDistributionType,
                          AKANTU_RANDOM_DISTRIBUTION_TYPES)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(RandomDistributionType,
                                AKANTU_RANDOM_DISTRIBUTION_TYPES)
AKANTU_CLASS_ENUM_INPUT_STREAM(RandomDistributionType,
                               AKANTU_RANDOM_DISTRIBUTION_TYPES)

/* -------------------------------------------------------------------------- */
/* Generator                                                                  */
/* -------------------------------------------------------------------------- */
template <typename T> class RandomGenerator {
  /* ------------------------------------------------------------------------ */
private:
  static long int _seed;                       // NOLINT
  static std::default_random_engine generator; // NOLINT
  /* ------------------------------------------------------------------------ */
public:
  inline T operator()() { return generator(); }

  /// function to print the contain of the class
  void printself(std::ostream & stream, int /* indent */) const {
    stream << "RandGenerator [seed=" << _seed << "]";
  }


  
  std::default_random_engine & getGenerator(); // Needed on apple clang
  /* ------------------------------------------------------------------------ */
public:
  static void seed(long int s) {
    _seed = s;
    generator.seed(_seed);
  }
  static long int seed() { return _seed; }

  static constexpr T min() { return std::default_random_engine::min(); }
  static constexpr T max() { return std::default_random_engine::max(); }
};

#if defined(__clang__)
template <typename T> long int RandomGenerator<T>::_seed; // NOLINT
template <typename T> std::default_random_engine RandomGenerator<T>::generator;
#endif

/* -------------------------------------------------------------------------- */
/* Some Helper                                                                */
/* -------------------------------------------------------------------------- */
template <typename T, class Distribution> struct RandomDistributionTypeHelper {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_not_defined;
};

template <>
struct RandomDistributionTypeHelper<Real,
                                    std::uniform_real_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_uniform;
};

template <>
struct RandomDistributionTypeHelper<Real, std::exponential_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_exponential;
};
template <>
struct RandomDistributionTypeHelper<Real, std::gamma_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_gamma;
};
template <>
struct RandomDistributionTypeHelper<Real, std::weibull_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_weibull;
};
template <>
struct RandomDistributionTypeHelper<Real,
                                    std::extreme_value_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_extreme_value;
};
template <>
struct RandomDistributionTypeHelper<Real, std::normal_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_normal;
};
template <>
struct RandomDistributionTypeHelper<Real, std::lognormal_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_lognormal;
};
template <>
struct RandomDistributionTypeHelper<Real, std::chi_squared_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_chi_squared;
};
template <>
struct RandomDistributionTypeHelper<Real, std::cauchy_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_cauchy;
};
template <>
struct RandomDistributionTypeHelper<Real, std::fisher_f_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_fisher_f;
};
template <>
struct RandomDistributionTypeHelper<Real, std::student_t_distribution<Real>> {
  static constexpr RandomDistributionType value =
      RandomDistributionType::_student_t;
};

/* -------------------------------------------------------------------------- */
template <class T> class RandomDistribution {
public:
  virtual ~RandomDistribution() = default;

  RandomDistribution() = default;
  RandomDistribution(const RandomDistribution & other) = default;
  RandomDistribution(RandomDistribution && other) noexcept = default;
  RandomDistribution & operator=(const RandomDistribution & other) = default;
  RandomDistribution &
  operator=(RandomDistribution && other) noexcept = default;

  virtual T operator()(RandomGenerator<Idx> & gen) = 0;
  virtual std::unique_ptr<RandomDistribution<T>> make_unique() const = 0;
  virtual void printself(std::ostream & stream, int = 0) const = 0;
};

template <class T, class Distribution>
class RandomDistributionProxy : public RandomDistribution<T> {
public:
  explicit RandomDistributionProxy(Distribution dist)
      : distribution(std::move(dist)) {}

  T operator()(RandomGenerator<Idx> & gen) override {
    return this->distribution(gen.getGenerator());
  }

  std::unique_ptr<RandomDistribution<T>> make_unique() const override {
    return std::make_unique<RandomDistributionProxy<T, Distribution>>(
        distribution);
  }

  void printself(std::ostream & stream, int /* indent */ = 0) const override {
    stream << std::to_string(
        RandomDistributionTypeHelper<T, Distribution>::value);
    stream << " [ " << distribution << " ]";
  }

private:
  Distribution distribution;
};

/* -------------------------------------------------------------------------- */
/* RandomParameter                                                            */
/* -------------------------------------------------------------------------- */
template <typename T> class RandomParameter {
public:
  template <class Distribution>
  explicit RandomParameter(T base_value, Distribution dist)
      : base_value(base_value),
        type(RandomDistributionType(
            RandomDistributionTypeHelper<T, Distribution>::value)),
        distribution_proxy(
            std::make_unique<RandomDistributionProxy<T, Distribution>>(
                std::move(dist))) {}

  explicit RandomParameter(T base_value)
      : base_value(base_value),
        type(RandomDistributionType(
            RandomDistributionTypeHelper<
                T, std::uniform_real_distribution<T>>::value)),
        distribution_proxy(
            std::make_unique<
                RandomDistributionProxy<T, std::uniform_real_distribution<T>>>(
                std::uniform_real_distribution<T>(0., 0.))) {}

  RandomParameter(const RandomParameter & other)
      : base_value(other.base_value), type(other.type),
        distribution_proxy(other.distribution_proxy->make_unique()) {}

  RandomParameter & operator=(const RandomParameter & other) {
    distribution_proxy = other.distribution_proxy->make_unique();
    base_value = other.base_value;
    type = other.type;
    return *this;
  }

  RandomParameter(RandomParameter && other) noexcept = default;
  RandomParameter & operator=(RandomParameter && other) noexcept = default;

  virtual ~RandomParameter() = default;

  inline void setBaseValue(const T & value) { this->base_value = value; }
  inline T getBaseValue() const { return this->base_value; }

  template <template <typename> class Generator, class iterator>
  void setValues(iterator it, iterator end) {
    RandomGenerator<Idx> gen;
    for (; it != end; ++it) {
      *it = this->base_value + (*distribution_proxy)(gen);
    }
  }

  virtual void printself(std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {
    stream << base_value;
    stream << " + " << *distribution_proxy;
  }

private:
  /// Value with no random variations
  T base_value;

  /// Random distribution type
  RandomDistributionType type;

  /// Proxy to store a std random distribution
  std::unique_ptr<RandomDistribution<T>> distribution_proxy;
};

/* -------------------------------------------------------------------------- */
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 RandomDistribution<T> & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 RandomParameter<T> & _this) {
  _this.printself(stream);
  return stream;
}

template <typename T>
auto make_random_parameter(
    T base,
    RandomDistributionType distribution = RandomDistributionType::_uniform,
    T a = T(), T b = T()) {
  switch (distribution) {
  case RandomDistributionType::_not_defined:
    return RandomParameter<T>(base, std::uniform_real_distribution<T>(0., 0.));
  case RandomDistributionType::_uniform:
    return RandomParameter<T>(base, std::uniform_real_distribution<T>(a, b));
  case RandomDistributionType::_exponential:
    return RandomParameter<T>(base, std::exponential_distribution<T>(a));
  case RandomDistributionType::_gamma:
    return RandomParameter<T>(base, std::gamma_distribution<T>(a, b));
  case RandomDistributionType::_weibull:
    return RandomParameter<T>(base, std::weibull_distribution<T>(b, a));
  case RandomDistributionType::_extreme_value:
    return RandomParameter<T>(base, std::extreme_value_distribution<T>(a, b));
  case RandomDistributionType::_normal:
    return RandomParameter<T>(base, std::normal_distribution<T>(a, b));
  case RandomDistributionType::_lognormal:
    return RandomParameter<T>(base, std::lognormal_distribution<T>(a, b));
  case RandomDistributionType::_chi_squared:
    return RandomParameter<T>(base, std::chi_squared_distribution<T>(a));
  case RandomDistributionType::_cauchy:
    return RandomParameter<T>(base, std::cauchy_distribution<T>(a, b));
  case RandomDistributionType::_fisher_f:
    return RandomParameter<T>(base, std::fisher_f_distribution<T>(a, b));
  case RandomDistributionType::_student_t:
    return RandomParameter<T>(base, std::student_t_distribution<T>(a));
  }

  return RandomParameter<T>(base, std::uniform_real_distribution<T>(0., 0.));
}



  template<typename T>
  std::default_random_engine & RandomGenerator<T>::getGenerator() {  return generator; } // Needed on apple clang
} // namespace akantu

#endif /* AKANTU_AKA_RANDOM_GENERATOR_HH_ */
