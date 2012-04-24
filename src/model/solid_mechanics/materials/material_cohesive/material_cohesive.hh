/**
 * @file   material_cohesive.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Feb  7 17:50:23 2012
 *
 * @brief  Specialization of the material class for cohesive elements
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
#include "material.hh"
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_COHESIVE_HH__
#define __AKANTU_MATERIAL_COHESIVE_HH__

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MaterialCohesive : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  typedef FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> > MyFEMCohesiveType;
public:

  MaterialCohesive(Model & model, const ID & id = "");
  virtual ~MaterialCohesive();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// read properties
  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// compute the residual for this material
  virtual void updateResidual(Vector<Real> & current_position,
			      GhostType ghost_type = _not_ghost);


  /// compute the stable time step for an element of size h
  virtual Real getStableTimeStep(__attribute__((unused)) Real h,
				 __attribute__((unused)) const Element & element = ElementNull) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  // /// add an element to the local mesh filter
  // __aka_inline__ void addElement(const ElementType & type,
  // 			 UInt element,
  // 			 const GhostType & ghost_type);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:

  /// constitutive law
  virtual void computeStress(__attribute__((unused)) ElementType el_type,
			     __attribute__((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute the tangent stiffness matrix
  virtual void computeTangentStiffness(__attribute__((unused)) const ElementType & el_type,
				       __attribute__((unused)) Vector<Real> & tangent_matrix,
				       __attribute__((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }


  void computeNormal(const Vector<Real> & position,
		     Vector<Real> & normal,
		     ElementType type,
		     GhostType ghost_type);

  void computeOpening(const Vector<Real> & displacement,
		      Vector<Real> & normal,
		      ElementType type,
		      GhostType ghost_type);

  template<ElementType type>
  void computeNormal(const Vector<Real> & position,
		     Vector<Real> & normal,
		     GhostType ghost_type);


  /// assemble residual
  void assembleResidual(GhostType ghost_type = _not_ghost);

  /// compute tractions (including normals and openings)
  void computeTraction(GhostType ghost_type = _not_ghost);

  /// constitutive law
  virtual void computeTraction(const Vector<Real> & normal,
			       ElementType el_type,
			       GhostType ghost_type = _not_ghost) = 0;

  /// compute reversible and total energies by element
  void computeEnergies();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// compute reversible energy
  Real getReversibleEnergy();

  /// compute dissipated energy
  Real getDissipatedEnergy();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// reversible energy by quadrature point
  ByElementTypeReal reversible_energy;

  /// total energy by quadrature point
  ByElementTypeReal total_energy;

  /// traction in all elements and quadrature points (previous time step)
  ByElementTypeReal tractions_old;

  /// opening in all elements and quadrature points (previous time step)
  ByElementTypeReal opening_old;

protected:

  /// traction in all elements and quadrature points
  ByElementTypeReal tractions;

  /// opening in all elements and quadrature points
  ByElementTypeReal opening;

  /// Link to the cohesive fem object in the model
  MyFEMCohesiveType * fem_cohesive;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_cohesive_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialCohesive & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_HH__ */
