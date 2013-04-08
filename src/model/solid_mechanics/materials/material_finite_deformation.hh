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

    /// initialize the material computed parameter
    virtual void initMaterial();

    /// compute the residual for this material
    //virtual void updateResidual(GhostType ghost_type = _not_ghost);

    /// interpolate   stress  on   given   positions  for   each  element   (empty
    /// implemantation to avoid the generic call to be done on cohesive elements)

    virtual void interpolateStress(__attribute__((unused)) const ElementType type,
            __attribute__((unused)) Array<Real> & result) {
    };

    /// assemble residual
    virtual void assembleResidual(GhostType ghost_type);
    
    virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

    /// assemble stiffness
    virtual void assembleStiffnessMatrix(GhostType ghost_type);

    template<UInt dim>
    void assembleStiffnessMatrixNL(const ElementType & type,
            GhostType ghost_type);
    
    template<UInt dim>
    void assembleStiffnessMatrixL2(const ElementType & type,
            GhostType ghost_type);

    
protected:
    
    /// stresses arrays ordered by element types
    ByElementTypeReal delta_stress;

    /// strains arrays ordered by element types
    ByElementTypeReal delta_strain;

    /// stresses arrays ordered by element types
    ByElementTypeReal stress_at_t;
    
    virtual void UpdateStressesAtT(GhostType ghost_type);

    
    inline UInt getCauchyStressMatrixSize(UInt spatial_dimension) const;
    inline UInt getCauchyStressArraySize(UInt spatial_dimension) const;
    
    template<UInt dim>
    inline void GreenStrain(const Matrix<Real> & F, Matrix<Real> & E);
    
    template<UInt dim>
    inline void AlmansiStrain(const Matrix<Real> & F, Matrix<Real> & E);
    
    template<UInt dim>
    inline void GreenStrain(const Matrix<Real> & d_u, const Matrix<Real> & u, Matrix<Real> & E);

    template<UInt dim>
    inline void transferBMatrixToBNL(const Matrix<Real> & B,
            Matrix<Real> & Bvoigt,
            UInt nb_nodes_per_element) const;
    
    template<UInt dim>
    inline void transferBMatrixToBL2(const Matrix<Real> & B, const Matrix<Real> & grad_u,
            Matrix<Real> & Bvoigt,
            UInt nb_nodes_per_element) const;
    
    template<UInt dim>
    inline void transferBMatrixToBL1(const Matrix<Real> & B,
            Matrix<Real> & Bvoigt,
            UInt nb_nodes_per_element) const;

    template<UInt dim>
    inline void SetCauchyStressMatrix(const Matrix<Real> & S_t, Matrix<Real> & Stress_matrix);
    
    void SetCauchyStressMatrix(const ElementType & type, Array<Real> & stress_matrix, GhostType ghost_type = _not_ghost);
    
    template<UInt dim>
    inline void SetCauchyStressArray(const Matrix<Real> & S_t, Matrix<Real> & Stress_vect);
    
    template<UInt dim>
    inline void SetCauchyStrainArray(const Matrix<Real> & E_t, Matrix<Real> & Strain_vect);
    
    void SetCauchyStressArray(const ElementType & type, Array<Real> & stress_vect, GhostType ghost_type = _not_ghost);


};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_finite_deformation_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_HH__ */
        