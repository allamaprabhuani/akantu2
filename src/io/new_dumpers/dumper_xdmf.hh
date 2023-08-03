/**
 * @file   dumper_xdmf.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Oct 03 2017
 *
 * @brief XDMF file writer
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
#include "dumper_hdf5.hh"
#include "dumper_file_base.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DUMPER_XDMF_HH__
#define __AKANTU_DUMPER_XDMF_HH__

namespace akantu {

class DumperXdmf : public DumperHDF5 {
public:
  explicit DumperXdmf(dumper::SupportBase & support);

protected:
  void dumpInternal() override;

  std::unique_ptr<dumper::FileBase> xdmf;
};

} // namespace akantu

#endif /* __AKANTU_DUMPER_XDMF_HH__ */
