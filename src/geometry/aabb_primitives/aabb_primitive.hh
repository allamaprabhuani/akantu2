/**
 * @file   aabb_primitive.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Fri Mar 13 2015
 *
 * @brief  Macro classe (primitive) for AABB CGAL algos
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_AABB_PRIMITIVE_HH__
#define __AKANTU_AABB_PRIMITIVE_HH__

#include "aka_common.hh"
#include "triangle.hh"
#include "tetrahedron.hh"

__BEGIN_AKANTU__

/**
 * This macro defines a class that is used in the CGAL AABB tree algorithm.
 * All the `typedef`s and methods are required by the AABB module.
 *
 * The member variables are
 *  - the id of the element associated to the primitive
 *  - the geometric primitive of the element
 *
 *  @param name the name of the primitive type
 *  @param kernel the name of the kernel used
 */
#define AKANTU_AABB_CLASS(name, kernel) \
  class name##_primitive {      \
    typedef std::list< name<kernel> >::iterator Iterator; \
                                                          \
  public:                                                 \
    typedef UInt Id;                                      \
    typedef kernel::Point_3 Point;                        \
    typedef kernel::name##_3 Datum;                       \
                                                          \
  public:                                                 \
    name##_primitive() : meshId(0), primitive() {}        \
    name##_primitive(Iterator it) : meshId(it->id()), primitive(*it) {} \
                                                                        \
  public:                                                               \
    const Datum & datum() const { return primitive; }                   \
    const Point & reference_point() const { return primitive.vertex(0); } \
    const Id & id() const { return meshId; }                              \
                                                                          \
  protected:                                                              \
    Id meshId;                                                            \
    name<kernel> primitive;                                               \
                                                                          \
  }

// If the primitive is supported by CGAL::intersection() then the 
// implementation process is really easy with this macro
AKANTU_AABB_CLASS(Triangle, CGAL::Cartesian<Real>);

#undef AKANTU_AABB_CLASS

__END_AKANTU__

#endif // __AKANTU_AABB_PRIMITIVE_HH__
