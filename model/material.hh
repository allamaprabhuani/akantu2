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
#include "aka_common.hh"
#include "aka_memory.hh"
#include "fem.hh"
#include "solid_mechanics_model.hh"


/* -------------------------------------------------------------------------- */
// namespace akantu {
//   class SolidMechanicsModel;
// };

__BEGIN_AKANTU__

class Material : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Material(SolidMechanicsModel & model, const MaterialID & id = "");
  virtual ~Material();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// read properties
  virtual void setParam(const std::string & key, const std::string & value) = 0;

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// constitutive law
  virtual void constitutiveLaw(ElementType el_type) = 0;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(ID, id, const MaterialID &);
  AKANTU_GET_MACRO(Rho, rho, Real);

  inline void setPotentialEnergyFlagOn();
  inline void setPotentialEnergyFlagOff();

  inline const Vector<Real> & getPotentialEnergy(ElementType type) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  friend class SolidMechanicsModel;

  /// id of the material
  MaterialID id;

  /// The model to witch the material belong
  SolidMechanicsModel * model;

  /// density : rho
  Real rho;

  /// list of element handled by the material
  ByConnectivityTypeUInt element_filter;

  /// has to compute potential energy or not
  bool potential_energy_flag;

  /// is the vector for potential energy initialized
  bool potential_energy_vector;

  /// potential energy by element
  ByConnectivityTypeReal potential_energy;

  /// boolean to know if the material has been initialized
  bool is_init;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Material & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#include "materials/material_elastic.hh"

#endif /* __AKANTU_MATERIAL_HH__ */

