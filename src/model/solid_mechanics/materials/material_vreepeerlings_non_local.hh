/**
 * @file   material_vreepeerlings_non_local.hh
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Feb 24 16:01:10 2012
 *
 * @brief  
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_vreepeerlings.hh"
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH__

__BEGIN_AKANTU__

/**
 * Material VreePeerlings Non local
 *
 * parameters in the material files :
 */
template<UInt spatial_dimension>
class MaterialVreePeerlingsNonLocal : public MaterialVreePeerlings<spatial_dimension>,
				      public MaterialNonLocal<spatial_dimension, BaseWeightFunction> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialNonLocal<spatial_dimension, BaseWeightFunction> MaterialNonLocalParent;

  MaterialVreePeerlingsNonLocal(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialVreePeerlingsNonLocal() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law
  virtual void computeNonLocalStress(GhostType ghost_type = _not_ghost);
  virtual void computeNonLocalStress(Vector<Real> & equi_strainnl,
				     ElementType el_type,
				     GhostType ghost_type = _not_ghost);


  /// Compute the tangent stiffness matrix for implicit for a given type
  void computeTangentStiffness(__attribute__ ((unused)) const ElementType & type,
			       __attribute__ ((unused)) Vector<double> & tangent_matrix,
			       __attribute__ ((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

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

  ByElementTypeReal equi_strain;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_vreepeerlings_non_local_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH__ */
