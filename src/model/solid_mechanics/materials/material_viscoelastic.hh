/**
 * @file   material_viscoelastic.hh
 * @author Vlad Yastrebov <vladislav.yastrebov@epfl.ch> 
 * @date   Thu Feb 7 2012 
 *
 * @brief  Material Visco-elastic, based on Standard Solid rheological model, see [] J.C. Simo, T.J.R. Hughes, "Computational Inelasticity", Springer (1998), see Sections 10.2 and 10.3
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
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VISCOELASTIC_HH__
#define __AKANTU_MATERIAL_VISCOELASTIC_HH__

__BEGIN_AKANTU__

/**
 * Material viscoelastic (caughey condition) isotropic
 *
 * parameters in the material files :
 *   - rho : density (default: 0)
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 *   - eta : viscosity
 *   - Ev  : stiffness of the viscous element
 *   - h   : internal variable of the integral history
 */
class MaterialViscoElastic : public MaterialElastic {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */  
public:

  MaterialViscoElastic(Model & model, const ID & id = "");

  virtual ~MaterialViscoElastic() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// constitutive law for a given quadrature point
  //__aka_inline__ void computeStress(Real * F, Real * sigma);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(HistoryIntegral, history_integral, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(StressDev,  stress_dev, Real);

  AKANTU_GET_MACRO(EV, Ev, const Real&); 
  AKANTU_SET_MACRO(EV, Ev, Real &); 

  AKANTU_GET_MACRO(Eta, eta, const Real&); 
  AKANTU_SET_MACRO(Eta, eta, Real &); 

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// stress due to viscosity
//  ByElementTypeReal stress_viscosity;

  /// stress due to elasticity
//  ByElementTypeReal stress_elastic;

  /// viscosity, viscous elastic modulus
  Real eta, Ev;

  /// history of deviatoric stress
  ByElementTypeReal stress_dev;

  /// Internal variable: history integral
  ByElementTypeReal history_integral;
};

/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_elastic_caughey_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/// standard output stream operator
/*
inline std::ostream & operator <<(std::ostream & stream, const MaterialViscoElastic & _this)
{
  _this.printself(stream);
  return stream;
}
*/
__END_AKANTU__

#endif /* __AKANTU_MATERIAL_VISCOELASTIC_HH__ */

