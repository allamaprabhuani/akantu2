/**
 * @file   aka_random_generator.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Nov  9 14:53:08 2012
 *
 * @brief  generic random generator
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
#include "aka_common.hh"
#include "aka_vector.hh"
#include "static_communicator.hh"

#ifndef __AKANTU_AKA_RANDOM_GENERATOR_HH__
#define __AKANTU_AKA_RANDOM_GENERATOR_HH__

__BEGIN_AKANTU__

template<typename T>
class RandomGenerator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  RandomGenerator(long int seed = 0) : seed(seed) {
    if(seed == 0) this->seed = time(NULL) * (StaticCommunicator::getStaticCommunicator().whoAmI() + 1);
    srand48(seed);
  }

  virtual ~RandomGenerator() {}


  virtual void generate(Array<T> & vect) = 0;
  virtual void setParams(std::string value) = 0;

  void setSeed(long int seed) {
    this->seed = seed;
    srand48(seed);
  }

  inline T rand();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

    stream << space << "seed=" << seed;
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  long int seed;
};

template<>
inline Real RandomGenerator<Real>::rand() {
  return drand48();
}

/// standard output stream operator
template<typename T>
inline std::ostream & operator <<(std::ostream & stream, const RandomGenerator<T>  & _this)
{
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
/* Uniform generator                                                          */
/* -------------------------------------------------------------------------- */
template<typename T>
class UniformRandomGenerator :public RandomGenerator<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  UniformRandomGenerator(long int seed = 0) : RandomGenerator<T>(seed) { };

  virtual ~UniformRandomGenerator() {};

  void generate(Array<T> & vect) {
    UInt n = vect.getSize();
    for (UInt i = 0; i < n; ++i)
      vect(i) = base + increment * this->rand();
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void setParams(std::string value) {
    std::stringstream sstr(value);
    Real top;
    sstr >> base;
    sstr >> top;
    increment = top - base;
  }

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    stream << "Uniform [ ";
    RandomGenerator<T>::printself(stream, indent);
    stream << ", min=" << base << ", max=" << base+increment << "]";
  };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  T base;
  T increment;
};


/* -------------------------------------------------------------------------- */
/* Weibull generator                                                          */
/* -------------------------------------------------------------------------- */
template<typename T>
class WeibullRandomGenerator :public RandomGenerator<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  WeibullRandomGenerator(long int seed = 0) : RandomGenerator<T>(seed) { };

  virtual ~WeibullRandomGenerator() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void generate(Array<T> & vect) {
    UInt n = vect.getSize();
    T e = T(1) / m;
    for (UInt i = 0; i < n; ++i) {
      T r = this->rand();
      vect(i) = minimum + lambda * std::pow(- std::log(1. - r), e);
    }
  }

  virtual void setParams(std::string value) {
    std::stringstream sstr(value);
    sstr >> m;
    sstr >> lambda;
    sstr >> minimum;
  }

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    stream << "Weibull [ ";
    RandomGenerator<T>::printself(stream, indent);
    stream << ", scale=" << m
	   << ", shape=" << lambda
	   << ", minimum=" << minimum << "]";
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_SET_MACRO(Shape, m, Real);
  AKANTU_SET_MACRO(Scale, lambda, Real);
  AKANTU_SET_MACRO(Minimum, minimum, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// shape parameter or Weibull modulus
  T m;
  /// scale parameter
  T lambda;
  /// minimum value
  T minimum;
};


__END_AKANTU__

#endif /* __AKANTU_AKA_RANDOM_GENERATOR_HH__ */
