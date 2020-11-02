/**
 * @file   material_cohesive_linear.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_SEQUENTIAL_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_SEQUENTIAL_HH__

namespace akantu {

/**
 * Cohesive material linear damage for extrinsic case
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default:
 * 0)
 *   - G_cI      : fracture energy for mode I (default: 0)
 *   - G_cII     : fracture energy for mode II (default: 0)
 *   - penalty   : stiffness in compression to prevent penetration
 */
template <UInt spatial_dimension>
class MaterialCohesiveLinearSequential
    : public MaterialCohesiveLinear<spatial_dimension> {
  /* -----------------------------------------------------------------------*/
  /* Constructors/Destructors                                               */
  /* -----------------------------------------------------------------------*/
public:
  MaterialCohesiveLinearSequential(SolidMechanicsModel & model,
                                   const ID & id = "");

  /* -----------------------------------------------------------------------*/
  /* Methods                                                                */
  /* -----------------------------------------------------------------------*/
public:
  /// initialise material
  void initMaterial() override;

  /// check stress for cohesive elements' insertion
  void checkInsertion(bool check_only = false) override;

  /// identify a facet with the highest tensile stress complying with
  /// topological criterion
  std::tuple<UInt, Real, Real>
  findCriticalFacet(const ElementType & type_facet);

  /// duplicate single facet and insert cohesive if check_only!=false
  void insertSingleCohesiveElement(const ElementType & type_facet,
                                   UInt facet_nb, Real crack_nb,
                                   bool check_only);

protected:
  /* -----------------------------------------------------------------------*/
  /* Accessors                                                              */
  /* -----------------------------------------------------------------------*/

  /* -----------------------------------------------------------------------*/
  /* Class Members                                                          */
  /* -----------------------------------------------------------------------*/
protected:
  /// portion of elements with the stress above their threshold to be inserted
  Real insertion_threshold;

  /// normal stresses on facets
  FacetInternalField<Real> normal_stresses;

  /// normal stresses on facets
  FacetInternalField<Real> normal_tractions;

  /// normal stresses normalized by the stress limit
  FacetInternalField<Real> effective_stresses;

  /// crack number
  CohesiveInternalField<Real> crack_number;
};

/* -------------------------------------------------------------------------*/
/* inline functions                                                         */
/* -------------------------------------------------------------------------*/

} // namespace akantu

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_SEQUENTIAL_HH__ */
