/**
 * @file   material_finite_deformation.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Wed Feb 22 16:31:20 2012
 *
 * @brief  Specialization of the material class for finite deformation
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
#include "fem_template.hh"
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_FINITE_DEFORMATION_HH__
#define __AKANTU_MATERIAL_FINITE_DEFORMATION_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
    class SolidMechanicsPlasticityModel;
}

__BEGIN_AKANTU__

class MaterialFiniteDeformation : public Material {
    /* ------------------------------------------------------------------------ */
    /* Constructors/Destructors                                                 */
    /* ------------------------------------------------------------------------ */
public:

    MaterialFiniteDeformation(SolidMechanicsModel& model, const ID & id = "");
    virtual ~MaterialFiniteDeformation();

    /* ------------------------------------------------------------------------ */
    /* Methods                                                                  */
    /* ------------------------------------------------------------------------ */
public:

    /// compute the residual for this material
    //virtual void updateResidual(GhostType ghost_type = _not_ghost);

    /// interpolate   stress  on   given   positions  for   each  element   (empty
    /// implemantation to avoid the generic call to be done on cohesive elements)

  //    virtual void interpolateStress(__attribute__((unused)) const ElementType type,
  //        __attribute__((unused)) Array<Real> & result) {
  // };


protected:

  /*   template<UInt dim>
    inline void GreenStrain(const Matrix<Real> & F, Matrix<Real> & E);

    template<UInt dim>
    inline void deltaGreenStrain(const Matrix<Real> & d_u, const Matrix<Real> & u, Matrix<Real> & E);*/

  /*template<UInt dim>
    inline void transferBMatrixToBL1(const Matrix<Real> & B,
            Matrix<Real> & Bvoigt,
            UInt nb_nodes_per_element) const;*/

  //    void SetCauchyStressArray(const ElementType & type, Array<Real> & stress_vect, GhostType ghost_type = _not_ghost);


};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_finite_deformation_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_HH__ */
