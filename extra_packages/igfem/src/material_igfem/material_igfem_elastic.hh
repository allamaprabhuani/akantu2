/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_igfem.hh"
#include "plane_stress_toolbox.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_IGFEM_ELASTIC_HH_
#define AKANTU_MATERIAL_IGFEM_ELASTIC_HH_

namespace akantu {

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template <Int spatial_dimension>
class MaterialIGFEMElastic
    : public PlaneStressToolbox<spatial_dimension, MaterialIGFEM> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef PlaneStressToolbox<spatial_dimension, MaterialIGFEM> Parent;

public:
  MaterialIGFEMElastic(SolidMechanicsModel & model, const ID & id = "");
  MaterialIGFEMElastic(SolidMechanicsModel & model, Int dim, const Mesh & mesh,
                       FEEngine & fe_engine, const ID & id = "");

  virtual ~MaterialIGFEMElastic() {}

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initMaterial();
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  virtual void computeTangentModuli(ElementType el_type,
                                    Array<Real> & tangent_matrix,
                                    GhostType ghost_type = _not_ghost);

  /// compute the elastic potential energy
  virtual void computePotentialEnergy(ElementType el_type,
                                      GhostType ghost_type = _not_ghost);

  virtual void
  computePotentialEnergyByElement(ElementType type, UInt index,
                                  Vector<Real> & epot_on_quad_points);

  void updateElasticInternals(const Array<Element> & element_list);

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
                                  Matrix<Real> & sigma,
                                  __attribute__((unused)) const Real lambda,
                                  const Real mu) const;

  /// compute the tangent stiffness matrix for an element
  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent,
                                         __attribute__((unused))
                                         const Real lambda,
                                         const Real mu) const;

  static inline void computePotentialEnergyOnQuad(const Matrix<Real> & grad_u,
                                                  const Matrix<Real> & sigma,
                                                  Real & epot);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  virtual void onElementsAdded(const Array<Element> & element_list,
                               const NewElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// First Lamé coefficient
  IGFEMInternalField<Real> lambda;

  /// Second Lamé coefficient (shear modulus)
  IGFEMInternalField<Real> mu;

  /// Bulk modulus
  IGFEMInternalField<Real> kpa;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_igfem_elastic_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_IGFEM_ELASTIC_HH_ */
