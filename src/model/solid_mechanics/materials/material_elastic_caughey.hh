/**
 * @file   material_elastic_caughey.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed May  4 15:16:59 2011
 *
 * @brief  Material isotropic viscoelastic (according to the Caughey condition)
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "material.hh"
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ELASTIC_CAUGHEY_HH__
#define __AKANTU_MATERIAL_ELASTIC_CAUGHEY_HH__

__BEGIN_AKANTU__

/**
 * Material viscoelastic (caughey condition) isotropic
 *
 * parameters in the material files :
 *   - rho : density (default: 0)
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 *   - alpha : viscous ratio
 */
template<UInt spatial_dimension>
class MaterialElasticCaughey : public MaterialElastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialElasticCaughey(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialElasticCaughey() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the potential energy for all elements
  virtual void computePotentialEnergy(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// constitutive law for a given quadrature point
  //inline void computeStress(Real * F, Real * sigma);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(StressViscosity, stress_viscosity, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(StressElastic,   stress_elastic, Real);

  AKANTU_GET_MACRO(Alpha, alpha, Real);
  AKANTU_SET_MACRO(Alpha, alpha, Real);

  virtual Real getParam(const ID & param) const;
  virtual void setParam(const ID & param, Real value);


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// stress due to viscosity
  ByElementTypeReal stress_viscosity;

  /// stress due to elasticity
  ByElementTypeReal stress_elastic;

  /// viscous ratio
  Real alpha;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_elastic_caughey_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_ELASTIC_CAUGHEY_HH__ */
