/**
 * @file   heat_transfer_interface_model.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Thu Apr 13 2023
 * @date last modification: Thu Apr 13 2023
 *
 * @brief  Model of Heat Transfer expanded to heat diffusion along interfaces
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "heat_transfer_model.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_HEAT_TRANSFER_INTERFACE_MODEL_HH_
#define AKANTU_HEAT_TRANSFER_INTERFACE_MODEL_HH_

namespace akantu {

class HeatTransferInterfaceModel : public HeatTransferModel {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MyFEEngineCohesiveType =
      FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive>;
  using Parent = HeatTransferModel;

  HeatTransferInterfaceModel(Mesh & mesh, UInt dim = _all_dimensions,
                             const ID & id = "heat_transfer_interface_model",
                             std::shared_ptr<DOFManager> dof_manager = nullptr);

  ~HeatTransferInterfaceModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// read one material file to instantiate all the materials
  void readMaterials() override;

  /// initialize the model
  void initModel() override;

  void predictor() override;
  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// compute and get the stable time step
  Real getStableTimeStep() override;

  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

public:
  // /// calculate the lumped capacity vector for heat transfer problem
  // void assembleCapacityLumped() override;

public:
  /// compute the internal heat flux
  void assembleInternalHeatRate() override;

  /// assemble the conductivity matrix
  void assembleConductivityMatrix() override;

  // /// compute normals on cohesive elements
  // void computeNormal(const Array<Real> & position, Array<Real> & normal,
  //                    ElementType type, GhostType ghost_type = _not_ghost);

  // /// compute basis alliged with cohesive elements
  // void computeBasis(const Array<Real> & position, Array<Real> & basis,
  //                   ElementType type, GhostType ghost_type = _not_ghost);

  /// assemble the capacity matrix
  void assembleCapacity() override;

  /// compute the capacity on quadrature points
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

  /// get FEEngine for cohesive elements
  const ShapeLagrange<_ek_cohesive> & getShapeFunctionsCohesive();

public:
  // /// asign material properties to physical groups
  // void assignPropertyToPhysicalGroup(const std::string & property_name,
  //                                    const std::string & group_name,
  //                                    Real value);
  // /// asign material properties to physical groups
  // void assignPropertyToPhysicalGroup(const std::string & property_name,
  //                                    const std::string & group_name,
  //                                    Matrix<Real> cond_matrix);

private:
  /// compute the integrated longitudinal conductivity matrix (or scalar in
  /// 2D) for each quadrature point
  void computeKLongOnQuadPoints(GhostType ghost_type);

  /// compute the transversal conductivity scalar devided by opening
  void computeKTransOnQuadPoints(GhostType ghost_type);

  void computeGradAndDeltaT(GhostType ghost_type);

  /// compute internal heat rate along the crack and assemble it into internal
  /// forces
  void computeLongHeatRate(GhostType ghost_type);

  /// compute internal heat rate perpendicular to the crack and assemble it into
  /// internal forces
  void computeTransHeatRate(GhostType ghost_type);

  /// compute heat rate contributions due to the opening/closure rate
  /// (dw/dt)
  void computeInertialHeatRate(GhostType ghost_type);

  /// compute an averaging operator of size (nb_nodes_per_itp_type x
  /// nb_nodes_per_cohesive_type)
  Matrix<Real> getAveragingOperator(const ElementType & type,
                                    const UInt & nb_degree_of_freedom = 1);

  /// compute an differencing operator of size (nb_nodes_per_itp_type x
  /// nb_nodes_per_cohesive_type)
  Matrix<Real> getDifferencingOperator(const ElementType & type,
                                       const UInt & nb_degree_of_freedom = 1);

  // /// compute the thermal energy
  // Real computeThermalEnergyByNode();

protected:
  /// calculate the lumped capacity vector for heat transfer problem (w
  /// ghost type)
  void assembleCapacityLumped(GhostType ghost_type) override {
    AKANTU_TO_IMPLEMENT();
  };
  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumpers::Field>
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       UInt spatial_dimension, ElementKind kind) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(TransversalConductivity, transversal_conductivity, Real);
  AKANTU_GET_MACRO(LongitudinalConductivity, longitudinal_conductivity, Real);
  AKANTU_GET_MACRO(CapacityInCrack, capacity_in_crack, Real);
  AKANTU_GET_MACRO(DensityInCrack, density_in_crack, Real);
  AKANTU_GET_MACRO(DefaultOpening, default_opening, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TransversalConductivityOnQpoints,
                                         k_perp_over_w, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(LongitudinalConductivityOnQpoints,
                                         k_long_w, Real);
  /// get the aperture on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Opening, opening_on_qpoints, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Opening, opening_on_qpoints, Real);

public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  Real transversal_conductivity;
  Real longitudinal_conductivity;
  Real default_opening;
  Real capacity_in_crack;
  Real density_in_crack;

  /// crack opening interpolated on quadrature points
  ElementTypeMapArray<Real> opening_on_qpoints;

  /// opening rate at quadrature points
  ElementTypeMapArray<Real> opening_rate;

  /// longitudinal conductivity tensor (or a scalar in 2D) multiplied by opening
  /// on quadrature points
  ElementTypeMapArray<Real> k_long_w;

  /// transversal conductivity scalar devided by opening on quadrature points
  ElementTypeMapArray<Real> k_perp_over_w;

  /// @brief boolean enabling opening-rate term into internal heat rate
  bool use_opening_rate{false};

  std::unordered_map<GhostType, UInt> long_conductivity_release{{_not_ghost, 0},
                                                                {_ghost, 0}};

  std::unordered_map<GhostType, UInt> perp_conductivity_release{{_not_ghost, 0},
                                                                {_ghost, 0}};

  UInt crack_conductivity_matrix_release{UInt(-1)};

  UInt opening_release{0};

  /// cohesive elements synchronizer
  std::unique_ptr<ElementSynchronizer> cohesive_synchronizer;
};

/* -------------------------------------------------------------------------- */
inline UInt
HeatTransferInterfaceModel::getNbData(const Array<Element> & elements,
                                      const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;
  Array<Element>::const_iterator<Element> it = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if (el.kind() == _ek_cohesive)
      std::cout << "Cohesive syncronized" << std::endl;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  size += Parent::getNbData(elements, tag);
  AKANTU_DEBUG_OUT();
  return size;
}
} // namespace akantu

#endif /* AKANTU_HEAT_TRANSFER_INTERFACE_MODEL_HH_ */
