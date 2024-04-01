/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "dumper_internal_types.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DUMPER_FUNCTORS_H_
#define AKANTU_DUMPER_FUNCTORS_H_

namespace akantu {
namespace dumper {
  /* ------------------------------------------------------------------------ */
  template <class ElementFunction> struct ElementTypeMapArrayFunctor {
    ElementTypeMapArrayFunctor(ElementType type, GhostType ghost_type,
                               const ElementFunction & function)
        : type(type), ghost_type(ghost_type), function(function) {}

    [[nodiscard]] Int getNbComponent(Int nb_component) const {
      if constexpr (details::has_getNbComponent_member<ElementFunction>) {
        return function.getNbComponent(nb_component, type, ghost_type);
      } else {
        return nb_component;
      }
    }

    [[nodiscard]] Int size(Int old_size) const {
      if constexpr (details::has_size_member<ElementFunction>) {
        return function.size(old_size, type, ghost_type);
      } else {
        return old_size;
      }
    }

    template <class Derived>
    [[nodiscard]] decltype(auto)
    operator()(const Eigen::MatrixBase<Derived> & in) const {
      return function(in, type, ghost_type);
    }

    ElementType type;
    GhostType ghost_type;
    const ElementFunction & function;
  };

  /* ------------------------------------------------------------------------ */
  struct NodeGroupsGlobalNodes {
    NodeGroupsGlobalNodes(const Mesh & mesh) : mesh(mesh) {}

    template <class Derived>
    Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & nodes) const {
      Vector<Idx> new_nodes(nodes.size());
      for (auto && [nn, on] : zip(new_nodes, nodes)) {
        nn = mesh.getNodeGlobalId(on);
      }
      return new_nodes;
    }

  private:
    const Mesh & mesh;
  };

  /* ------------------------------------------------------------------------ */
  struct ElementGroupGlobalElements {
    ElementGroupGlobalElements(const ElementGroup & group) : group(group) {}

    template <class Derived>
    Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & elements,
                           ElementType type, GhostType ghost_type) const {
      Vector<Idx> new_elements(elements.size());
      auto offset = group.getMesh().getElementsOffsets(type, ghost_type);
      for (auto && [ne, oe] : zip(new_elements, elements)) {
        ne = offset + oe;
      }
      return new_elements;
    }

  private:
    const ElementGroup & group;
  };

  /* ------------------------------------------------------------------------ */
  struct toVTKConnectivity {
    toVTKConnectivity(const Mesh & mesh) : mesh(mesh) {
      if (not write_reorder_initialized) {
        for (auto type : element_types) {
          Vector<Idx> reorder(Mesh::getNbNodesPerElement(type));
          for (auto && [i, n] : enumerate(reorder)) {
            n = i;
          }

          switch (type) {
          case _cohesive_3d_12:
            reorder << 0, 1, 2, 6, 7, 8, 3, 4, 5, 9, 10, 11;
            break;
          case _cohesive_2d_6:
            reorder << 0, 2, 1, 4, 5, 3;
            break;
          case _cohesive_2d_4:
            reorder << 0, 1, 3, 2;
            break;
          case _hexahedron_20:
            reorder[12] = 16;
            reorder[13] = 17;
            reorder[14] = 18;
            reorder[15] = 19;
            reorder[16] = 12;
            reorder[17] = 13;
            reorder[18] = 14;
            reorder[19] = 15;
            break;
          case _pentahedron_15:
            reorder << 0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11;
            break;
          default:
            break;
          }
          write_reorder[type] = reorder;
        }
        write_reorder_initialized = true;
      }
    }

    template <class Derived>
    Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & connectivity,
                           ElementType type, GhostType /*ghost_type*/) const {
      Vector<Idx> new_connectivity(connectivity.size());
      auto && reorder = write_reorder[type];
      // Rewriting to global connectivities
      for (auto && [nc, c] : zip(new_connectivity, reorder)) {
        nc = mesh.getNodeGlobalId(connectivity(c));
      }
      return new_connectivity;
    }

  private:
    const Mesh & mesh;
    static std::map<ElementType, Vector<Idx>> write_reorder;
    static bool write_reorder_initialized;
  };

  /* ------------------------------------------------------------------------ */
  struct ElementGroupConnectivityFunctor {
    ElementGroupConnectivityFunctor(const Mesh & mesh)
        : rewrite(mesh), connectivities(mesh.getConnectivities()) {}

    [[nodiscard]] Int getNbComponent(Int /*nb_component*/, ElementType type,
                                     GhostType ghost_type) const {
      return connectivities(type, ghost_type).getNbComponent();
    }

    template <class Derived>
    Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & element,
                           ElementType type, GhostType ghost_type) const {
      Element el{type, element[0], ghost_type};
      return rewrite(connectivities.get(el), type, ghost_type);
    }

  private:
    toVTKConnectivity rewrite;
    const ElementTypeMapArray<Idx> & connectivities;
  };

  /* ------------------------------------------------------------------------ */
  struct ElementGroupNodesFunctor {
    ElementGroupNodesFunctor(const Array<Real> & nodes) : nodes(nodes) {}

    [[nodiscard]] Int getNbComponent(Int /*nb_component*/
    ) const {
      return nodes.getNbComponent();
    }

    template <class Derived>
    Vector<Real> operator()(const Eigen::MatrixBase<Derived> & node) const {
      return make_view(nodes, nodes.getNbComponent()).begin()[node[0]];
    }

  private:
    const Array<Real> & nodes;
  };
} // namespace dumper
} // namespace akantu

#endif // AKANTU_DUMPER_FUNCTORS_H_
