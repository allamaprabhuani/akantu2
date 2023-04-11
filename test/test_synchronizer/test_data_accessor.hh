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
#include "data_accessor.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;
/* -------------------------------------------------------------------------- */

class TestAccessor : public DataAccessor<Element> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  inline TestAccessor(const Mesh & mesh,
                      const ElementTypeMapArray<Real> & barycenters);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Barycenter, barycenters, Real);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  inline Int getNbData(const Array<Element> & elements,
                       const SynchronizationTag & tag) const;
  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  const ElementTypeMapArray<Real> & barycenters;
  const Mesh & mesh;
};

/* -------------------------------------------------------------------------- */
/* TestSynchronizer implementation                                            */
/* -------------------------------------------------------------------------- */
inline TestAccessor::TestAccessor(const Mesh & mesh,
                                  const ElementTypeMapArray<Real> & barycenters)
    : barycenters(barycenters), mesh(mesh) {}

inline Int TestAccessor::getNbData(const Array<Element> & elements,
                                   const SynchronizationTag &) const {
  return mesh.getSpatialDimension() * sizeof(Real) * elements.size();
}

inline void TestAccessor::packData(CommunicationBuffer & buffer,
                                   const Array<Element> & elements,
                                   const SynchronizationTag &) const {
  for (const auto & element : elements) {
    auto && bary = this->barycenters.get(element);
    buffer << bary;
  }
}

inline void TestAccessor::unpackData(CommunicationBuffer & buffer,
                                     const Array<Element> & elements,
                                     const SynchronizationTag &) {
  for (const auto & element : elements) {
    auto && barycenter_loc = this->barycenters.get(element);
    Vector<Real> bary(barycenter_loc.size());
    buffer >> bary;

    auto dist = (barycenter_loc - bary).template lpNorm<Eigen::Infinity>();
    EXPECT_NEAR(0, dist, 1e-15);
  }
}
