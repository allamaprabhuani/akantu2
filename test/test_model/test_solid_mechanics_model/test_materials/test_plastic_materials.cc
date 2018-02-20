/* -------------------------------------------------------------------------- */
#include "material_linear_isotropic_hardening.hh"
#include "solid_mechanics_model.hh"
#include "test_material_fixtures.hh"
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

using types = ::testing::Types<Traits<MaterialLinearIsotropicHardening, 1>,
                               Traits<MaterialLinearIsotropicHardening, 2>,
                               Traits<MaterialLinearIsotropicHardening, 3>>;

/*****************************************************************/

namespace {

template <typename T>
class TestPlasticMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_CASE(TestPlasticMaterialFixture, types);

TYPED_TEST(TestPlasticMaterialFixture, PlasticComputeStress) {
  this->material->testComputeStress();
}
TYPED_TEST(TestPlasticMaterialFixture, PlasticEnergyDensity) {
  this->material->testEnergyDensity();
}
TYPED_TEST(TestPlasticMaterialFixture, PlasticComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}
TYPED_TEST(TestPlasticMaterialFixture, PlasticComputeCelerity) {
  this->material->testCelerity();
}
}

/*****************************************************************/
