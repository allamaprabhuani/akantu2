/**
 * @file   test_material_fixtures.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Nov 17 2017
 * @date last modification:  Fri Dec 13 2019
 *
 * @brief  Fixture for material tests
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <material_elastic.hh>
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <random>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
template <typename T> class FriendMaterial : public T {
public:
  ~FriendMaterial() = default;

  FriendMaterial(SolidMechanicsModel & model, const ID & id = "material")
      : T(model, id) {
    gen.seed(::testing::GTEST_FLAG(random_seed));
  }

  virtual void testComputeStress() { AKANTU_TO_IMPLEMENT(); };
  virtual void testComputeTangentModuli() { AKANTU_TO_IMPLEMENT(); };
  virtual void testEnergyDensity() { AKANTU_TO_IMPLEMENT(); };
  virtual void testCelerity() { AKANTU_TO_IMPLEMENT(); }

  virtual void setParams() {}
  virtual void SetUp() {
    this->setParams();
    this->initMaterial();
  }

  virtual void TearDown() {}

  inline Matrix<Real> getDeviatoricStrain(Real intensity);
  inline Matrix<Real> getHydrostaticStrain(Real intensity);
  inline Matrix<Real> getComposedStrain(Real intensity);

  inline Matrix<Real> reverseRotation(const Matrix<Real> & A,
                                      const Matrix<Real> & R) {
    Matrix<Real> res = R.transpose() * A * R;
    return res;
  };

  inline Matrix<Real> applyRotation(const Matrix<Real> & A,
                                    const Matrix<Real> & R) {
    Matrix<Real> res = R * A * R.transpose();
    return res;
  };

  inline Vector<Real> getRandomVector();
  inline Matrix<Real> getRandomRotation();

protected:
  std::mt19937 gen;
};

/* -------------------------------------------------------------------------- */
template <typename T>
Matrix<Real> FriendMaterial<T>::getDeviatoricStrain(Real intensity) {
  Matrix<Real> dev{{0., 1., 2.}, {1., 0., 3.}, {2., 3., 0.}};
  dev *= intensity;
  return dev;
}

/* -------------------------------------------------------------------------- */
template <typename T>
Matrix<Real> FriendMaterial<T>::getHydrostaticStrain(Real intensity) {
  Matrix<Real> dev{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  dev *= intensity;
  return dev;
}

/* -------------------------------------------------------------------------- */

template <typename T>
Matrix<Real> FriendMaterial<T>::getComposedStrain(Real intensity) {
  Matrix<Real> s = FriendMaterial<T>::getHydrostaticStrain(intensity) +
                   FriendMaterial<T>::getDeviatoricStrain(intensity);
  s *= intensity;
  return s;
}

/* -------------------------------------------------------------------------- */
template <typename T> Vector<Real> FriendMaterial<T>::getRandomVector() {
  auto dim = this->spatial_dimension;
  std::uniform_real_distribution<Real> dis;
  Vector<Real> vector(dim);
  vector.zero();
  while (vector.norm() < 1e-9) {
    for (auto s : arange(dim)) {
      vector(s) = dis(gen);
    }
  }
  return vector;
}

/* -------------------------------------------------------------------------- */
template <typename T> Matrix<Real> FriendMaterial<T>::getRandomRotation() {
  const auto dim = this->spatial_dimension;

  Matrix<Real> rotation(dim, dim);
  auto && v1_ = rotation(0);
  v1_ = getRandomVector().normalized();
  if (dim == 2) {
    Vector<Real, 3> v1{v1_(0), v1_(1), 0};
    Vector<Real, 3> v3{0, 0, 1};

    Vector<Real, 3> v2 = v3.cross(v1);

    rotation(1) = v2.block<2, 1>(0, 0);
  }

  if (dim == 3) {
    Vector<Real, 3> v2 = getRandomVector();
    Vector<Real, 3> v1 = v1_;
    v2 = (v2 - v2.dot(v1) * v1).normalized();

    Vector<Real, 3> v3 = v1.cross(v2);

    rotation(1) = v2;
    rotation(2) = v3;
  }

//#define debug_
#if defined(debug_)
  if (dim == 2)
    rotation = Matrix<Rea, 2, 2>{{1., 0.}, {0., 1.}};
  if (dim == 3)
    rotation = Matrix<Real, 3, 3>{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
#endif

  rotation.transposeInPlace();

  return rotation;
}

/* -------------------------------------------------------------------------- */
template <typename T, class Model>
class TestMaterialBaseFixture : public ::testing::Test {
public:
  using mat_class = typename T::mat_class;

  void SetUp() override {
    constexpr auto spatial_dimension = T::Dim;
    mesh = std::make_unique<Mesh>(spatial_dimension);
    model = std::make_unique<Model>(*mesh);
    material = std::make_unique<friend_class>(*model);
    material->SetUp();
  }

  void TearDown() override {
    material->TearDown();
    material.reset(nullptr);
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  using friend_class = FriendMaterial<mat_class>;

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<Model> model;
  std::unique_ptr<friend_class> material;
};

/* -------------------------------------------------------------------------- */
template <template <Int> class T, Int _Dim> struct Traits {
  using mat_class = T<_Dim>;
  static constexpr Int Dim = _Dim;
};

/* -------------------------------------------------------------------------- */
template <typename T>
using TestMaterialFixture = TestMaterialBaseFixture<T, SolidMechanicsModel>;

/* -------------------------------------------------------------------------- */
