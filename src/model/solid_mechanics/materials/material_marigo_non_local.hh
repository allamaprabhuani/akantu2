/**
 * @file   material_marigo_non_local.hh
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Aug 31 17:08:23 2011
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
#include "aka_common.hh"
#include "material_marigo.hh"
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_MARIGO_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_MARIGO_NON_LOCAL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

/**
 * Material Marigo
 *
 * parameters in the material files :
 */
template<UInt spatial_dimension, template <UInt> class WeightFunction = BaseWeightFunction>
class MaterialMarigoNonLocal : public MaterialMarigo<spatial_dimension>,
			       public MaterialNonLocal<spatial_dimension, WeightFunction> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialNonLocal<spatial_dimension, WeightFunction> MaterialNonLocalParent;

  MaterialMarigoNonLocal(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialMarigoNonLocal() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);


  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law
  virtual void computeNonLocalStress(GhostType ghost_type = _not_ghost);

  virtual void computeNonLocalStress(Vector<Real> & damage,
				     ElementType el_type,
				     GhostType ghost_type = _not_ghost);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

private:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Y, Y, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ByElementTypeReal Y;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_marigo_non_local_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_MARIGO_NON_LOCAL_HH__ */
