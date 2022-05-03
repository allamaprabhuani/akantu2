/**
 * @file   solid_mechanics_model_RVE.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jan 13 14:54:18 2016
 * @update Tue Feb 8  2022
 *
 * @brief  SMM for RVE computations in FE2 simulations
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
#ifndef __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__

/* -------------------------------------------------------------------------- */
#include "aka_grid_dynamic.hh"
#include "rve_tools.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class SolidMechanicsModelRVE : public SolidMechanicsModel, public RVETools {

  /* ----------------------------------------------------------------- */
  /* Constructors/Destructors                                          */
  /* ----------------------------------------------------------------- */

public:
  SolidMechanicsModelRVE(Mesh & mesh, UInt nb_expanding_elements,
                         bool use_RVE_mat_selector = true,
                         UInt dim = _all_dimensions,
                         const ID & id = "solid_mechanics_model",
                         std::shared_ptr<DOFManager> dof_manager = nullptr);

  virtual ~SolidMechanicsModelRVE();

  typedef VoigtHelper<2> voigt_h;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  void initFullImpl(const ModelOptions & option) override;

  /// initialize the materials
  void initMaterials() override;

public:
  /// advance expansion by applying eigen strain at certain material and apply
  /// homogenized properties
  void advanceExpansion(const Matrix<Real> & prestrain,
                        const ID & material_name = "gel");

  /// correct the rigid boundary movement and assemble internal forces
  void assembleInternalForces() override;
  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<UInt> & index,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(CornerNodes, corner_nodes, const Array<UInt> &);
  AKANTU_GET_MACRO(Volume, volume, Real);
  bool stiffnessChanged() { return this->stiffness_changed; };

private:
  /* ------------------------------------------------------------------------ */
  /* Members */
  /* ------------------------------------------------------------------------ */
  /// standard mat selector or user one
  bool use_RVE_mat_selector;

  /// the number of expanding finite elements inside the RVE
  UInt nb_expanding_elements;

  /// the number of gel pockets inside the RVE
  bool stiffness_changed;
};

inline void SolidMechanicsModelRVE::unpackData(CommunicationBuffer & buffer,
                                               const Array<UInt> & index,
                                               const SynchronizationTag & tag) {
  SolidMechanicsModel::unpackData(buffer, index, tag);
}
} // namespace akantu

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__ */
