#include "aka_bounding_box.hh"
#include "aka_vector.hh"

__BEGIN_AKANTU__


template <>
std::ostream& operator<< <1>(std::ostream& os, const BoundingBox<1>& bb) {
  
  os<<"Line["<<bb.min()<<","<<bb.max()<<"]";
  return os;
}


template <>
std::ostream& operator<< <2>(std::ostream& os, const BoundingBox<2>& bb) {
    
    os<<"Rectangle["<<bb.min()<<","<<bb.max()<<"]";
    return os;    
}

template <>
std::ostream& operator<< <3>(std::ostream& os, const BoundingBox<3>& bb) {
    
    os<<"Cuboid["<<bb.min()<<","<<bb.max()<<"]";
    return os;
}




/* -------------------------------------------------------------------------- */
template <int dim, class nodes_container>
BoundingBox<dim> createPointList(const nodes_container& nodes, const Array<Real>& coord) {
    
    //  AKANTU_DEBUG_ASSERT(nodes.getSize() != 0, "No nodes to create a bounding box with.");
    typedef typename nodes_container::const_iterator node_iterator;
    
    node_iterator it = nodes.begin();
    assert(it != nodes.end());
    
    BoundingBox<dim> bbox(Point<dim>(&coord(*it),coord.getNbComponent()));
    for (++it; it != nodes.end(); ++it) {
        Real * point_coord = &coord(*it);
        for (UInt d=0; d<coord.getNbComponent(); ++d) {
            ;
        }
        bbox += *it;
    }
    return bbox;
}

__END_AKANTU__
