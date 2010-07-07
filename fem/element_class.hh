/**
 * @file   element_class.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jun 16 18:48:25 2010
 *
 * @brief
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_HH__
#define __AKANTU_ELEMENT_CLASS_HH__

/* -------------------------------------------------------------------------- */
#include "common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template<ElementType type> class ElementClass {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ElementClass();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Accesors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(NbNodesPerElement, nb_nodes_per_element, unsigned int);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  unsigned int nb_nodes_per_element;

  unsigned int nb_quadrature_points;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

template<ElementType type> ElementClass<type>::ElementClass() {
  nb_nodes_per_element = 0;
  nb_quadrature_points = 0;
}


/* -------------------------------------------------------------------------- */

__END_AKANTU__


#endif /* __AKANTU_ELEMENT_CLASS_HH__ */
