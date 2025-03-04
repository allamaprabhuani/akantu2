/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dumper_iohelper.hh"
/* -------------------------------------------------------------------------- */
#ifndef AKANTU_DUMPER_TEXT_HH_
#define AKANTU_DUMPER_TEXT_HH_
/* -------------------------------------------------------------------------- */
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

namespace akantu {

class DumperText : public DumperIOHelper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DumperText(const std::string & basename = "dumper_text",
             iohelper::TextDumpMode mode = iohelper::_tdm_space,
             bool parallel = true);
  ~DumperText() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void
  registerMesh(const Mesh & mesh, Int spatial_dimension = _all_dimensions,
               GhostType ghost_type = _not_ghost,
               ElementKind element_kind = _ek_not_defined) override;

  void registerFilteredMesh(
      const Mesh & mesh, const ElementTypeMapArray<Idx> & elements_filter,
      const Array<Idx> & nodes_filter,
      Int spatial_dimension = _all_dimensions,
      GhostType ghost_type = _not_ghost,
      ElementKind element_kind = _ek_not_defined) override;

  void setBaseName(const std::string & basename) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  void setPrecision(Int prec);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // namespace akantu

#endif /* AKANTU_DUMPER_TEXT_HH_ */
