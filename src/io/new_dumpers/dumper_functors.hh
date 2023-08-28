#include "aka_array.hh"
#include "dumper_internal_types.hh"

#ifndef AKANTU_DUMPER_FUNCTORS_H_
#define AKANTU_DUMPER_FUNCTORS_H_

namespace akantu {
namespace dumper {
  /* ------------------------------------------------------------------------ */
  template <class ElementFunction> struct ElementTypeMapArrayFunctor {
    ElementTypeMapArrayFunctor(ElementType type,
                               const ElementFunction & function)
        : type(type), function(function) {}

    [[nodiscard]] Int getNbComponent(Int nb_component) const {
      if constexpr (details::has_getNbComponent_member<ElementFunction>) {
        return function.getNbComponent(nb_component, type);
      } else {
        return nb_component;
      }
    }

    [[nodiscard]] Int size(Int old_size) const {
      if constexpr (details::has_size_member<ElementFunction>) {
        return function.size(old_size, type);
      } else {
        return old_size;
      }
    }

    template <class Derived>
    [[nodiscard]] decltype(auto)
    operator()(const Eigen::MatrixBase<Derived> & in) const {
      return function(in, type);
    }

    ElementType type;
    const ElementFunction & function;
  };

  /* ------------------------------------------------------------------------ */
  struct toVTKConnectivity {
    template <class Derived>
    Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & connectivity,
                           ElementType /*type*/) const {
      return connectivity;
    }
  };

  /* ------------------------------------------------------------------------ */
  struct ElementGroupConnectivityFunctor {
    ElementGroupConnectivityFunctor(
        const ElementTypeMapArray<Idx> & connectivities)
        : connectivities(connectivities) {}

    [[nodiscard]] Int getNbComponent(Int /*nb_component*/,
                                     ElementType type) const {
      return connectivities(type).getNbComponent();
    }

    template <class Derived>
    Vector<Idx> operator()(const Eigen::MatrixBase<Derived> & element,
                           ElementType type) const {
      Element el{type, element[0], _not_ghost};
      return rewrite(connectivities.get(el), type);
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
