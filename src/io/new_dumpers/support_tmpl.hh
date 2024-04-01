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
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "node_group.hh"
/* -------------------------------------------------------------------------- */
#include <map>
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

      this->positions = make_field(mesh.getNodes(), *this);
      this->connectivities =
          make_field(mesh.getConnectivities(), *this, toVTKConnectivity(mesh));

      positions->addProperty("name", "position");
      connectivities->addProperty("name", "connectivities");

      this->addProperty("name", mesh.getID());
    }

    void updateOffsets() override { mesh.updateOffsets(); }

    [[nodiscard]] ElementTypesIteratorHelper
    elementTypes(GhostType ghost_type = _not_ghost) const override {
      return mesh.elementTypes(mesh.getSpatialDimension(), ghost_type,
                               _ek_not_defined);
    }

    [[nodiscard]] Release getRelease() const override {
      return mesh.getRelease();
    }

    [[nodiscard]] Int getNbNodes() const override { return mesh.getNbNodes(); }
    [[nodiscard]] Int getNbLocalNodes() const { return mesh.getNbLocalNodes(); }
    [[nodiscard]] Int getNbGlobalNodes() const override {
      return mesh.getNbGlobalNodes();
    }

    [[nodiscard]] Int
    getNbElements(const ElementType & type,
                  const GhostType & ghost_type = _not_ghost) const override {
      return mesh.getConnectivities()(type, ghost_type).size();
    }
    [[nodiscard]] Int getNbGlobalElements(
        const ElementType & type,
        const GhostType & ghost_type = _not_ghost) const override {
      return mesh.getNbGlobalElements(type, ghost_type);
    }

    [[nodiscard]] Int getElementsOffsets(
        const ElementType & type,
        const GhostType & ghost_type = _not_ghost) const override {
      return mesh.getElementsOffsets(type, ghost_type);
    }

    [[nodiscard]] bool isDistributed() const override {
      return mesh.isDistributed();
    }

    [[nodiscard]] const Communicator & getCommunicator() const override {
      return mesh.getCommunicator();
    }

    Mesh & getMesh() { return mesh; }
    const Mesh & getMesh() const { return mesh; }

  protected:
    Mesh & mesh;
  };

  /* ---------------------------------------------------------------------- */
  template <>
  class Support<ElementGroup> : public SupportBase, public SupportElements {
  public:
    explicit Support(ElementGroup & inner)
        : SupportBase(SupportType::_element_group), inner(inner) {

      auto && mesh = inner.getMesh();

      this->positions = make_field(inner.getNodeGroup().getNodes(), *this,
                                   ElementGroupNodesFunctor(mesh.getNodes()));
      this->connectivities = make_field(inner.getElements(), *this,
                                        ElementGroupConnectivityFunctor(mesh));
      positions->addProperty("name", "positions");
      connectivities->addProperty("name", "connectivities");

      nodes_list = make_field(inner.getNodeGroup().getNodes(), *this,
                              NodeGroupsGlobalNodes(mesh));

      elements = make_field(inner.getElements(), *this,
                            ElementGroupGlobalElements(inner));

      nodes_list->addProperty("name", "nodes");
      elements->addProperty("name", "elements");

      this->addProperty("name", inner.getName());
    }

    //    void updateOffsets() override { inner.getMesh().updateOffsets(); }

    [[nodiscard]] auto & getMesh() const { return inner.getMesh(); }

    [[nodiscard]] auto & getElements() const { return *elements; }
    [[nodiscard]] auto & getNodeIds() const { return *nodes_list; }

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

    [[nodiscard]] bool isDistributed() const override {
      return inner.getMesh().isDistributed();
    }

    [[nodiscard]] Int getNbNodes() const override {
      return inner.getNodeGroup().getNbNodes();
    }

    [[nodiscard]] Int getNbLocalNodes() const override {
      return inner.getNodeGroup().getNbLocalNodes();
    }

    [[nodiscard]] Int getNbGlobalNodes() const override {
      return inner.getNodeGroup().getNbGlobalNodes();
    }

    [[nodiscard]] Int
    getNbElements(const ElementType & type,
                  const GhostType & ghost_type = _not_ghost) const override {
      return inner.getElements()(type, ghost_type).size();
    }

    [[nodiscard]] Int getNbGlobalElements(
        const ElementType & type,
        const GhostType & ghost_type = _not_ghost) const override {
      return inner.getNbGlobalElements(type, ghost_type);
    }

    [[nodiscard]] Int getElementsOffsets(
        const ElementType & type,
        const GhostType & ghost_type = _not_ghost) const override {
      return inner.getElementsOffsets(type, ghost_type);
    }

    [[nodiscard]] Int getNodesOffsets() const {
      return inner.getNodeGroup().getOffsets();
    }

  protected:
    ElementGroup & inner;
    std::shared_ptr<FieldElementMapArrayTemplateBase<Idx>> elements;
    std::shared_ptr<FieldNodalArrayTemplateBase<Idx>> nodes_list;
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
