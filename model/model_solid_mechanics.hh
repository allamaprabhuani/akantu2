/**
 * @file   model_solid_mechanics.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 11:51:06 2010
 *
 * @brief  Model of Solid Mechanics
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __MYFEM_MODEL_SOLID_MECHANICS_HH__
#define __MYFEM_MODEL_SOLID_MECHANICS_HH__

__BEGIN_AKANTU__

class ModelSolidMechanics : Model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  ModelSolidMechanics();
  virtual ~ModelSolidMechanics();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "model_solid_mechanics_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const ModelSolidMechanics & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __MYFEM_MODEL_SOLID_MECHANICS_HH__ */
