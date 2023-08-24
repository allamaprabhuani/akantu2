/**
 * @file   dumper_xdmf.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Oct 03 2017
 *
 * @brief Dump data in xdmf format
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
#include "dumper_xdmf.hh"
#include "aka_types.hh"
#include "support.hh"
#include "xdmf_file.hh"
/* -------------------------------------------------------------------------- */
#include <filesystem>
/* -------------------------------------------------------------------------- */

namespace fs = std::filesystem;

namespace akantu {

struct XdmfConnectivityFunctor {
  template <class Derived>
  Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & connectivity,
                         ElementType type) const {
    Eigen::PermutationMatrix<Eigen::Dynamic> P(connectivity.size());
    P.setIdentity();

    auto & perm = P.indices();
    switch (type) {
    case _triangle_6:
      perm[3] = 5;
      perm[4] = 3;
      perm[5] = 4;
      break;
    case _tetrahedron_10:
      perm[4] = 9;
      perm[5] = 6;
      perm[6] = 8;
      perm[7] = 7;
      perm[8] = 5;
      perm[9] = 4;
      break;
    case _quadrangle_4:
      perm[2] = 3;
      perm[3] = 2;
      break;
    default:
      // nothing to change
      break;
    }

    return P.transpose() * connectivity;
  }
};

/* -------------------------------------------------------------------------- */
DumperXdmf::DumperXdmf(dumper::SupportBase & support) : DumperHDF5(support) {}

/* -------------------------------------------------------------------------- */
void DumperXdmf::dumpInternal() {
  if (time_activated) {
    support.addProperty("dump_count", count);
    support.addProperty("time", Real(count));
  }

  if (aka::is_of_type<dumper::Support<Mesh>>(support)) {
    auto add_xdmf_connectivities_for_support = [](auto && support) {
      if (not support.hasField("xdmf_connectivities")) {
        auto xdmf_connectivities =
            dumper::make_field(aka::as_type<dumper::SupportElements>(support)
                                   .getConnectivities()
                                   .getSharedPointer(),
                               support, XdmfConnectivityFunctor());
        xdmf_connectivities->addProperty("xdmf_ignore_field", true);
        support.addField("xdmf_connectivities", xdmf_connectivities);
      }
    };

    add_xdmf_connectivities_for_support(support);

    for (auto && [_, sub_support] : support.getSubSupports()) {
      add_xdmf_connectivities_for_support(*sub_support);
    }
  }

  DumperHDF5::dumpInternal();
  using dumper::XDMF::File;

  if (not xdmf) {
    auto path = fs::path(directory);
    path /= filename + ".xmf";
    xdmf = std::make_unique<File>(support, path);
  }

  xdmf->dump();
}

} // namespace akantu
