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

#ifndef __MYFEM_ELEMENT_CLASS_HH__
#define __MYFEM_ELEMENT_CLASS_HH__

/* -------------------------------------------------------------------------- */
#include "common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_MYFEM__

template<ElementType type> class ElementClass {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ElementClass();
  virtual ~ElementClass();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accesors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

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

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const ElementClass & _this)
{
  _this.printself(stream);
}


__END_MYFEM__


#endif /* __MYFEM_ELEMENT_CLASS_HH__ */
