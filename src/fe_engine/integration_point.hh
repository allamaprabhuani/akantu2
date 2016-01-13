/**
 * @file   integration_point.hh
 *
 * @author 
 *
 * @date creation: 
 * @date last modification: Mon Oct 19 2015
 *
 * @brief  definition of the class IntegrationPoint
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_QUADRATURE_POINT_H
#define AKANTU_QUADRATURE_POINT_H
/* -------------------------------------------------------------------------- */
#include "element.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */
class IntegrationPoint;
extern const IntegrationPoint IntegrationPointNull;
/* -------------------------------------------------------------------------- */

class IntegrationPoint : public Element {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:
  typedef Vector<Real> position_type;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  IntegrationPoint(const Element & element, UInt num_point = 0, UInt nb_quad_per_element = 0) :
    Element(element), num_point(num_point),
    global_num(element.element*nb_quad_per_element + num_point),
    position((Real *)NULL, 0) { };

  IntegrationPoint(ElementType type = _not_defined, UInt element = 0,
		  UInt num_point = 0, GhostType ghost_type = _not_ghost) :
    Element(type, element, ghost_type), num_point(num_point), global_num(0),
    position((Real *)NULL, 0) { };

  IntegrationPoint(UInt element, UInt num_point,
		  UInt global_num,
		  const position_type & position,
		  ElementType type,
		  GhostType ghost_type = _not_ghost) :
    Element(type, element, ghost_type), num_point(num_point), global_num(global_num),
    position((Real *)NULL, 0) { this->position.shallowCopy(position); };

  IntegrationPoint(const IntegrationPoint & quad) :
    Element(quad), num_point(quad.num_point), global_num(quad.global_num), position((Real *) NULL, 0) {
    position.shallowCopy(quad.position);
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  inline bool operator==(const IntegrationPoint & quad) const {
    return Element::operator==(quad) && this->num_point == quad.num_point;
  }

  inline bool operator!=(const IntegrationPoint & quad) const {
    return ((element != quad.element)
            || (type != quad.type)
            || (ghost_type != quad.ghost_type)
            || (kind != quad.kind)
	    || (num_point != quad.num_point)
	    || (global_num != quad.global_num));
  }

  bool operator<(const IntegrationPoint& rhs) const {
    bool res = Element::operator<(rhs) || (Element::operator==(rhs) && this->num_point < rhs.num_point);
    return res;
  }

  inline IntegrationPoint & operator=(const IntegrationPoint & q) {
    if(this != &q) {
      element    = q.element;
      type       = q.type;
      ghost_type = q.ghost_type;
      num_point  = q.num_point;
      global_num = q.global_num;
      position.shallowCopy(q.position);
    }

    return *this;
  }
  
  /// get the position of the integration point
  AKANTU_GET_MACRO(Position, position, const position_type &);

  /// set the position of the integration point
  void setPosition(const position_type & position) {
    this->position.shallowCopy(position);
  }

  /// deep copy of the position of the integration point
  void copyPosition(const position_type & position) {
    this->position.deepCopy(position);
  }

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
    stream << space << "IntegrationPoint [";
    Element::printself(stream, 0);
    stream << ", " << num_point << "(" << global_num << ")" << "]";
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

public:
  /// number of quadrature point in the element
  UInt num_point;
  /// global number of the quadrature point
  UInt global_num;
  // TODO might be temporary: however this class should be tought maybe...
  std::string material_id;
private:

  /// position of the quadrature point
  position_type position;
};

__END_AKANTU__


#endif /* AKANTU_QUADRATURE_POINT_H */
