/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_DUMPER_FIELD_HH_
#define AKANTU_DUMPER_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_iohelper.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumpers {
  /* ------------------------------------------------------------------------ */
  class FieldComputeProxy;
  class FieldComputeBaseInterface;
  class ComputeFunctorInterface;
  class HomogenizerProxy;
  /* ------------------------------------------------------------------------ */

  /// Field interface
  class Field : public std::enable_shared_from_this<Field> {
    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors */
    /* ---------------------------------------------------------------------- */
  public:
    Field() = default;
    virtual ~Field() = default;

    /* ---------------------------------------------------------------------- */
    /* Methods */
    /* ---------------------------------------------------------------------- */
  public:
    /// register this to the provided dumper
    virtual void registerToDumper(const std::string & id,
                                  iohelper::Dumper & dumper) = 0;

    /// set the number of data per item (used for elements fields at the moment)
    virtual void setNbData(Int /*nb_data*/) { AKANTU_TO_IMPLEMENT(); };

    /// set the number of data per elem (used for elements fields at the moment)
    virtual void
    setNbDataPerElem([[gnu::unused]] const ElementTypeMap<Int> & nb_data) {
      AKANTU_TO_IMPLEMENT();
    };

    /// set the number of data per elem (used for elements fields at the moment)
    virtual void setNbDataPerElem([[gnu::unused]] Int nb_data) {
      AKANTU_TO_IMPLEMENT();
    };

    /// get the number of components of the hosted field
    virtual ElementTypeMap<Int>
    getNbComponents(Int /*dim*/ = _all_dimensions,
                    GhostType /*ghost_type*/ = _not_ghost,
                    ElementKind /*kind*/ = _ek_not_defined) {
      AKANTU_TO_IMPLEMENT();
    };

    /// for connection to a FieldCompute
    inline virtual std::shared_ptr<Field>
    connect(FieldComputeProxy & /*proxy*/) {
      AKANTU_TO_IMPLEMENT();
    };

    /// for connection to a FieldCompute
    inline virtual std::unique_ptr<ComputeFunctorInterface>
    connect(HomogenizerProxy & /*proxy*/) {
      AKANTU_TO_IMPLEMENT();
    };

    /// check if the same quantity of data for all element types
    virtual void checkHomogeneity() = 0;

    /// return the dumper name
    std::string getGroupName() { return group_name; };

    /// return the id of the field
    std::string getID() { return field_id; };

    /* ---------------------------------------------------------------------- */
    /* Accessors */
    /* ---------------------------------------------------------------------- */
  public:
    /// return the flag to know if the field is homogeneous/contiguous
    virtual bool isHomogeneous() { return homogeneous; }

    /* ---------------------------------------------------------------------- */
    /* Class Members */
    /* ---------------------------------------------------------------------- */
  protected:
    /// the flag to know if it is homogeneous
    bool homogeneous{false};

    /// the name of the group it was associated to
    std::string group_name;

    /// the name of the dumper it was associated to
    std::string field_id;
  };

  /* ------------------------------------------------------------------------ */

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_FIELD_HH_ */
