/**
 * @file   material_elastic.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 29 15:00:59 2010
 *
 * @brief  Material isotropic elastic
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ELASTIC_HH__
#define __AKANTU_MATERIAL_ELASTIC_HH__

__BEGIN_AKANTU__

class MaterialElastic : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialElastic(SolidMechanicsModel & model, const MaterialID & id = "");

  virtual ~MaterialElastic() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  void setParam(const std::string & key, const std::string & value) {  };

  /// constitutive law for all element of a type
  void constitutiveLaw(ElementType el_type);

  /// constitutive law for a given quadrature point
  inline void constitutiveLaw(Real * F, Real * sigma, Real * epot);

  /// compute the celerity of wave in the material
  inline Real celerity();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// the young modulus
  Real E;

  /// Poisson coefficient
  Real nu;

  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// Bulk modulus
  Real kpa;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_elastic_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialElastic & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_ELASTIC_HH__ */
