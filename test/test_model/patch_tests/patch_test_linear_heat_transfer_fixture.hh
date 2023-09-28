/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "heat_transfer_model.hh"
/* -------------------------------------------------------------------------- */
#include "patch_test_linear_fixture.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PATCH_TEST_LINEAR_HEAT_TRANSFER_FIXTURE_HH_
#define AKANTU_PATCH_TEST_LINEAR_HEAT_TRANSFER_FIXTURE_HH_

/* -------------------------------------------------------------------------- */
template <typename type>
class TestPatchTestHTMLinear
    : public TestPatchTestLinear<type, HeatTransferModel> {
  using parent = TestPatchTestLinear<type, HeatTransferModel>;

public:
  void applyBC() override {
    parent::applyBC();
    auto & temperature = this->model->getTemperature();
    this->applyBConDOFs(temperature);
  }

  void initModel(const AnalysisMethod & method,
                 const std::string & material_file) override {
    TestPatchTestLinear<type, HeatTransferModel>::initModel(method,
                                                            material_file);
    if (method != _static) {
      this->model->setTimeStep(0.5 * this->model->getStableTimeStep());
    }
  }

  void checkAll() {
    auto & temperature = this->model->getTemperature();
    auto && cl = this->model->getConstitutiveLaw(0);
    Matrix<Real> C = cl.get("conductivity");
    this->checkDOFs(temperature);

    auto && grad_T = cl.template getArray<Real>("∇u", this->type);
    this->checkGradient(grad_T, temperature);

    auto && K_grad_T = cl.template getArray<Real>("D∇u", this->type);

    this->checkResults(
        [&](const Matrix<Real> & grad_T) { return C * grad_T.transpose(); },
        K_grad_T, temperature);
  }
};

using htm_types = gtest_list_t<TestElementTypes>;

TYPED_TEST_SUITE(TestPatchTestHTMLinear, htm_types, );

#endif /* AKANTU_PATCH_TEST_LINEAR_HEAT_TRANSFER_FIXTURE_HH_ */
