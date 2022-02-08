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
#include "asr_tools.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
//#include <unordered_set>
/* -------------------------------------------------------------------------- */

namespace akantu {

class SolidMechanicsModelRVE : public SolidMechanicsModel, public ASRTools {

  /* ----------------------------------------------------------------- */
  /* Constructors/Destructors                                          */
  /* ----------------------------------------------------------------- */

public:
  SolidMechanicsModelRVE(Mesh & mesh, bool use_RVE_mat_selector = true,
                         UInt nb_gel_pockets = 400, UInt dim = _all_dimensions,
                         const ID & id = "solid_mechanics_model",
                         const MemoryID & memory_id = 0);

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
  /// advance the reactions -> grow gel and apply homogenized properties
  void advanceASR(const Matrix<Real> & prestrain);

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

  /// the number of gel pockets inside the RVE
  UInt nb_gel_pockets;

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
