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
#include "support.hh" // NOLINT(pp_including_mainfile_in_preamble)
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

  /* ------------------------------------------------------------------------ */
  template <> class Support<Mesh> : public SupportBase, public SupportElements {
  public:
    using ElementTypesIteratorHelper =
        SupportElements::ElementTypesIteratorHelper;

    explicit Support(Mesh & mesh)
        : SupportBase(SupportType::_mesh), mesh(mesh) {

      this->nodes = make_field(mesh.getNodes(), *this);
      this->connectivities =
          make_field(mesh.getConnectivities(), *this, toVTKConnectivity());

      nodes->addProperty("name", "position");
      connectivities->addProperty("name", "connectivities");

      this->addProperty("name", mesh.getID());
    }

    [[nodiscard]] ElementTypesIteratorHelper
    elementTypes(GhostType ghost_type = _not_ghost) const override {
      return mesh.elementTypes(mesh.getSpatialDimension(), ghost_type,
                               _ek_not_defined);
    }

    [[nodiscard]] Release getRelease() const override {
      return mesh.getRelease();
    }

    [[nodiscard]] Int getNbNodes() const override { return mesh.getNbNodes(); }
    [[nodiscard]] Int
    getNbElements(const ElementType & type,
                  const GhostType & ghost_type = _not_ghost) const override {
      return mesh.getConnectivities()(type, ghost_type).size();
    }

  protected:
    Mesh & mesh;
  };

  /* ------------------------------------------------------------------------ */
  template <>
  class Support<ElementGroup> : public SupportBase, public SupportElements {
  public:
    explicit Support(ElementGroup & inner)
        : SupportBase(SupportType::_element_group), inner(inner),
          elements(make_field(inner.getElements(), *this)),
          nodes_list(make_field(inner.getNodeGroup().getNodes(), *this)) {

      this->nodes =
          make_field(inner.getNodeGroup().getNodes(), *this,
                     ElementGroupNodesFunctor(inner.getMesh().getNodes()));
      this->connectivities = make_field(
          inner.getElements(), *this,
          ElementGroupConnectivityFunctor(inner.getMesh().getConnectivities()));

      nodes->addProperty("name", "nodes");
      elements->addProperty("name", "elements");
      connectivities->addProperty("name", "xdmf_connectivities");

      this->addProperty("name", inner.getName());
    }

    [[nodiscard]] auto & getElements() const { return *elements; }

    [[nodiscard]] ElementTypesIteratorHelper
    elementTypes(GhostType ghost_type = _not_ghost) const override {
      return inner.elementTypes(inner.getMesh().getSpatialDimension(),
                                ghost_type, _ek_not_defined);
    }

    [[nodiscard]] Release getRelease() const override {
      return inner.getMesh().getRelease();
    }

    [[nodiscard]] bool contains(ElementType type,
                                GhostType ghost_type = _not_ghost) const {
      auto && types = elementTypes(ghost_type);
      return (std::find(types.begin(), types.end(), type) != types.end());
    }

    [[nodiscard]] Int getNbNodes() const override {
      return inner.getNodeGroup().getNodes().size();
    }
    [[nodiscard]] Int
    getNbElements(const ElementType & type,
                  const GhostType & ghost_type = _not_ghost) const override {
      return inner.getElements()(type, ghost_type).size();
    }

  protected:
    ElementGroup & inner;
    std::shared_ptr<FieldElementMapArrayTemplateBase<Idx>> elements;
    std::shared_ptr<FieldArrayTemplateBase<Idx>> nodes_list;
  };

  /* ------------------------------------------------------------------------
   */
  template <typename T> auto & SupportBase::addSubSupport(T & type) {
    auto && ptr = make_support(type);
    auto & support = *ptr;
    ptr->parent = this;
    sub_supports_[support.getName()] = std::move(ptr);
    return support;
  }

  /* ------------------------------------------------------------------------
   */
  const auto & SupportBase::getSubSupports() const { return sub_supports_; }
  const auto & SupportBase::getParentSupport() const { return *parent; }

  /* ------------------------------------------------------------------------
   */
  inline void SupportBase::addField(const ID & id,
                                    const std::shared_ptr<FieldBase> & field,
                                    FieldUsageType usage) {
    field->addProperty("name", id);
    field->addProperty("field_usage", usage);
    fields_.emplace(id, field);
  }

  /* ------------------------------------------------------------------------
   */
  template <class Cont, std::enable_if_t<not std::is_convertible_v<
                            std::decay_t<Cont>, std::shared_ptr<FieldBase>>> *>
  inline void SupportBase::addField(const ID & id, Cont && data,
                                    FieldUsageType usage) {
    auto && field = make_field(std::forward<Cont>(data), *this);
    this->addField(id, std::move(field), usage);
  }

  /* ------------------------------------------------------------------------
   */
  template <class Cont, class Func,
            std::enable_if_t<not std::is_convertible_v<
                std::decay_t<Cont>, std::shared_ptr<FieldBase>>> *>
  inline void SupportBase::addField(const ID & id, Cont && data, Func && func,
                                    FieldUsageType usage) {
    auto && field =
        make_field(std::forward<Cont>(data), *this, std::forward<Func>(func));
    this->addField(id, std::move(field), usage);
  }

} // namespace dumper
} // namespace akantu

#endif /* __AKANTU_SUPPORT_TMPL_HH__ */
