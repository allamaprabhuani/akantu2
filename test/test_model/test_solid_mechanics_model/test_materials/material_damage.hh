/**
 * @file   material_elastic.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 29 15:00:59 2010
 *
 * @brief  Material isotropic elastic
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DAMAGE_HH__
#define __AKANTU_MATERIAL_DAMAGE_HH__

__BEGIN_AKANTU__

class MaterialDamage : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialDamage(SolidMechanicsModel & model, const MaterialID & id = "");

  virtual ~MaterialDamage() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  void setParam(const std::string & key, const std::string & value,
		const MaterialID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law for a given quadrature point
  inline void computeStress(Real * F, Real * sigma,Real & damage);

  /// compute the potential energy for all elements
  void computePotentialEnergy(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the potential energy for on element
  inline void computePotentialEnergy(Real * F, Real * sigma, Real * epot);

  /// compute the celerity of wave in the material
  inline Real celerity();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the stable time step
  inline Real getStableTimeStep(Real h);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Damage, damage, const Vector<Real> &);
private:

  /// the young modulus
  Real E;

  /// Poisson coefficient
  Real nu;

  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// resistance to damage
  Real Yd;

  /// damage threshold
  Real Sd;

  /// Bulk modulus
  Real kpa;

  /// damage internal variable
  ByElementTypeReal damage;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_damage_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialDamage & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_DAMAGE_HH__ */
