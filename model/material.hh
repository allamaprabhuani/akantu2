/**
 * @file   material.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jul 23 09:06:29 2010
 *
 * @brief  Mother class for all materials
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_HH__
#define __AKANTU_MATERIAL_HH__

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "memory.hh"
//#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  class SolidMechanicsModel;
};

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/*  MaterialBase                                                              */
/* -------------------------------------------------------------------------- */

class MaterialBase : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialBase(SolidMechanicsModel & model, const MaterialID & id = "");
  virtual ~MaterialBase() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// read properties
  void SetParam(const std::string & key, const std::string & value) { };

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const { };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(Rho, rho, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// id of the material
  MaterialID id;

  /// The model to witch the material belong
  SolidMechanicsModel * model;

  /// density : rho
  Real rho;
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*  Material                                                                  */
/* -------------------------------------------------------------------------- */

template<MaterialType type>
class Material : public MaterialBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Material(SolidMechanicsModel & model, const MaterialID & id = "");
  virtual ~Material() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
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

#include "material_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialBase & _this)
{
  _this.printself(stream);
  return stream;
}

/// standard output stream operator
template<MaterialType type>
inline std::ostream & operator <<(std::ostream & stream, const Material<type> & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_HH__ */
