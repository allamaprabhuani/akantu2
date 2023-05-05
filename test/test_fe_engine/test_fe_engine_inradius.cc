/**
 * Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <fstream>
/* -------------------------------------------------------------------------- */
#include "test_fe_engine_fixture.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

TYPED_TEST(TestFEMFixture, Inradius) {
  const auto & connectivities = this->mesh->getConnectivity(this->type);
  auto nodes = this->mesh->getNodes().cbegin(this->dim);

  std::ifstream fin;
  fin.open(std::to_string(this->type) + ".csv");
  if (not fin.is_open()) {
    throw std::runtime_error("Could not open file " +
                             std::to_string(this->type) + ".csv");
  }

  int id;
  std::vector<Real> inradiuses(connectivities.size());
  Real inradius;

  std::string line;
  while (std::getline(fin, line)) {
    // Create a stringstream of the current line
    std::stringstream ss(line);
    ss >> id;
    // If the next token is a comma, ignore it and move on
    if (ss.peek() == ',') {
      ss.ignore();
    }
    ss >> inradius;

    inradiuses[id] = inradius;
  }

  Matrix<Real> X(this->dim, connectivities.getNbComponent());
  Element element{this->type, 0, _not_ghost};
  for (auto && conn :
       make_view(connectivities, connectivities.getNbComponent())) {
    for (auto s : arange(conn.size())) {
      X(s) = nodes[conn(s)];
    }

    auto inradius_coor = this->fem->getElementInradius(X, this->type);
    auto inradius_elem = this->fem->getElementInradius(element);
    auto inradius_file = inradiuses[element.element];

    // Files was written in simple precision...
    EXPECT_NEAR(inradius_elem, inradius_file, 1e-6);
    EXPECT_NEAR(inradius_coor, inradius_file, 1e-6);

    ++element.element;
  }
}
