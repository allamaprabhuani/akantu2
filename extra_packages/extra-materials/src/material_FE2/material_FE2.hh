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
#include "material.hh"
#include "material_thermal.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_FE_2_HH_
#define AKANTU_MATERIAL_FE_2_HH_

namespace akantu {
class SolidMechanicsModelRVE;
}

namespace akantu {

/* -------------------------------------------------------------------------- */
/// /!\ This material works ONLY for meshes with a single element type!!!!!
/* -------------------------------------------------------------------------- */

/**
 * MaterialFE2
 *
 * parameters in the material files :
 *   - mesh_file
 */
template <UInt DIM> class MaterialFE2 : public MaterialThermal<DIM> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef MaterialThermal<DIM> Parent;

public:
  MaterialFE2(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialFE2();

  typedef VoigtHelper<DIM> voigt_h;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initMaterial();

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

  /// advance alkali-silica reaction
  void advanceASR(const Matrix<Real> & prestrain);

private:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Underlying RVE at each integration point
  std::vector<std::unique_ptr<SolidMechanicsModelRVE>> RVEs;

  /// Meshes for all RVEs
  std::vector<std::unique_ptr<Mesh>> meshes;

  /// the element type of the associated mesh (this material handles only one
  /// type!!)
  ElementType el_type;

  /// the name of RVE mesh file
  ID mesh_file;

  /// Elastic stiffness tensor at each Gauss point (in voigt notation)
  InternalField<Real> C;

  /// number of gel pockets in each underlying RVE
  UInt nb_gel_pockets;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_FE2_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_FE_2_HH_ */
