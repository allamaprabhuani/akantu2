/**
 * @file   support_tmpl.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Fri Oct 06 2017
 *
 * @brief A Documented file.
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
#include "support.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_field.hh"
#include "element_group.hh"
#include "mesh_accessor.hh"
#include "node_group.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SUPPORT_TMPL_HH__
#define __AKANTU_SUPPORT_TMPL_HH__

namespace akantu {
namespace dumper {

  template <typename T> auto & SupportBase::addSubSupport(const T & type) {
    auto && ptr = make_support(type);
    auto & support = *ptr;
    ptr->parent = this;
    sub_supports_[support.getName()] = std::move(ptr);
    return support;
  }

  const auto & SupportBase::getSubSupports() const { return sub_supports_; }
  const auto & SupportBase::getParentSupport() const { return *parent; }

  void SupportBase::addField(const ID & id,
                             std::unique_ptr<FieldBase> && field) {
    field->addProperty("name", id);
    fields_.emplace(id, std::move(field));
  }

  /* ------------------------------------------------------------------------ */
  template <> class Support<Mesh> : public SupportBase, public SupportElements {
  public:
    using ElementTypesIteratorHelper =
        SupportElements::ElementTypesIteratorHelper;

    explicit Support(Mesh & inner)
        : SupportBase(SupportType::_mesh), inner(inner),
          nodes(make_field(MeshAccessor(inner).getNodes(), *this)),
          connectivities(make_field(MeshAccessor(inner).getConnectivities(),
                                    *this, ConnectivityFunctor())) {

      nodes->addProperty("name", "position");
      connectivities->addProperty("name", "connectivities");

      this->addProperty("name", inner.getID());
      this->addProperty("release", this->getRelease());
    }

    [[nodiscard]] FieldNodeArray<Real> & getNodes() const { return *nodes; }
    [[nodiscard]] FieldConnectivity & getConnectivities() const {
      return *connectivities;
    };

    [[nodiscard]] ElementTypesIteratorHelper
    elementTypes(GhostType ghost_type = _not_ghost) const override {
      return inner.elementTypes(inner.getSpatialDimension(), ghost_type,
                                _ek_not_defined);
    }

    [[nodiscard]] Int getRelease() { return inner.getRelease(); }

  protected:
    Mesh & inner;
    std::unique_ptr<FieldNodeArray<Real>> nodes;
    std::unique_ptr<FieldConnectivity> connectivities;
  };

  /* ------------------------------------------------------------------------ */
  template <>
  class Support<ElementGroup> : public SupportBase, public SupportElements {
  public:
    explicit Support(ElementGroup & inner)
        : SupportBase(SupportType::_element_group), inner(inner),
          nodes(make_field(inner.getNodeGroup().getNodes(), *this)),
          elements(make_field(inner.getElements(), *this)) {

      nodes->addProperty("name", "nodes");
      elements->addProperty("name", "elements");

      this->addProperty("name", inner.getName());
      this->addProperty("release", this->getRelease());
    }

    [[nodiscard]] auto & getNodes() const { return *nodes; }
    [[nodiscard]] auto & getElements() const { return *elements; };

    [[nodiscard]] ElementTypesIteratorHelper
    elementTypes(GhostType ghost_type = _not_ghost) const override {
      return inner.elementTypes(inner.getMesh().getSpatialDimension(),
                                ghost_type, _ek_not_defined);
    }

    [[nodiscard]] Int getRelease() { return inner.getMesh().getRelease(); }

  protected:
    ElementGroup & inner;
    std::unique_ptr<FieldNodeArray<Idx>> nodes;
    std::unique_ptr<FieldElementMapArray<Idx>> elements;
  };

} // namespace dumper
} // namespace akantu

#endif /* __AKANTU_SUPPORT_TMPL_HH__ */
