/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_HEAT_TRANSFER_MODEL_HH_
#define AKANTU_HEAT_TRANSFER_MODEL_HH_

namespace akantu {
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

namespace akantu {

class HeatTransferModel : public Model,
                          public DataAccessor<Element>,
                          public DataAccessor<Idx> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  HeatTransferModel(Mesh & mesh, Int spatial_dimension = _all_dimensions,
                    const ID & id = "heat_transfer_model",
                    std::shared_ptr<DOFManager> dof_manager = nullptr);

  ~HeatTransferModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// read one material file to instantiate all the materials
  void readMaterials();

  /// allocate all vectors
  void initSolver(TimeStepSolverType time_step_solver_type,
                  NonLinearSolverType non_linear_solver_type) override;

  /// initialize the model
  void initModel() override;

  void predictor() override;

  /// compute the heat flux
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) const override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID & matrix_id) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;
  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// compute and get the stable time step
  Real getStableTimeStep();

  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  // temporary protection to prevent bad usage: should check for bug
protected:
  /// compute the internal heat flux \todo Need code review: currently not
  /// public method
  void assembleInternalHeatRate();

public:
  /// calculate the lumped capacity vector for heat transfer problem
  void assembleCapacityLumped();

public:
  /// assemble the conductivity matrix
  void assembleConductivityMatrix();

  /// assemble the conductivity matrix
  void assembleCapacity();

  /// compute the capacity on quadrature points
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

private:
  /// calculate the lumped capacity vector for heat transfer problem (w
  /// ghost type)
  void assembleCapacityLumped(GhostType ghost_type);

  /// compute the conductivity tensor for each quadrature point in an array
  void computeConductivityOnQuadPoints(GhostType ghost_type);

  /// compute vector \f[k \grad T\f] for each quadrature point
  void computeKgradT(GhostType ghost_type);

  /// compute the thermal energy
  Real computeThermalEnergyByNode();

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline Int getNbData(const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;
  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  inline Int getNbData(const Array<Idx> & indexes,
                       const SynchronizationTag & tag) const override;
  inline void packData(CommunicationBuffer & buffer, const Array<Idx> & indexes,
                       const SynchronizationTag & tag) const override;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Idx> & indexes,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumpers::Field>
  createNodalFieldReal(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createNodalFieldBool(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       Int spatial_dimension, ElementKind kind) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_AUTO(Density, density);
  AKANTU_GET_MACRO_AUTO(Capacity, capacity);
  /// get the dimension of the system space
  AKANTU_GET_MACRO_AUTO(SpatialDimension, spatial_dimension);
  /// get the current value of the time step
  AKANTU_GET_MACRO_AUTO(TimeStep, time_step);
  /// get the assembled heat flux
  AKANTU_GET_MACRO_DEREF_PTR(InternalHeatRate, internal_heat_rate);
  /// get the boundary vector
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(BlockedDOFs, blocked_dofs);
  /// get the external heat rate vector
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ExternalHeatRate, external_heat_rate);
  /// get the temperature gradient
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureGradient,
                                         temperature_gradient, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConductivityOnQpoints,
                                         conductivity_on_qpoints, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureOnQpoints,
                                         temperature_on_qpoints, Real);
  /// internal variables
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(KgradT, k_gradt_on_qpoints, Real);
  /// get the temperature
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Temperature, temperature);
  /// get the temperature derivative
  AKANTU_GET_MACRO_DEREF_PTR(TemperatureRate, temperature_rate);

  /// get the energy denominated by thermal
  Real getEnergy(const std::string & energy_id, ElementType type, Idx index);
  /// get the energy denominated by thermal
  Real getEnergy(const std::string & energy_id);

  /// get the thermal energy for a given element
  Real getThermalEnergy(ElementType type, Idx index);
  /// get the thermal energy for a given element
  Real getThermalEnergy();

protected:
  /* ------------------------------------------------------------------------ */
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ----------------------------------------------------------------------- */
  template <class iterator, class const_iterator>
  void getThermalEnergy(iterator Eth, const_iterator T_it,
                        const_iterator T_end) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// time step
  Real time_step;

  /// temperatures array
  std::unique_ptr<Array<Real>> temperature;

  /// temperatures derivatives array
  std::unique_ptr<Array<Real>> temperature_rate;

  /// increment array (@f$\delta \dot T@f$ or @f$\delta T@f$)
  std::unique_ptr<Array<Real>> increment;

  /// the density
  Real density;

  /// the speed of the changing temperature
  ElementTypeMapArray<Real> temperature_gradient;

  /// temperature field on quadrature points
  ElementTypeMapArray<Real> temperature_on_qpoints;

  /// conductivity tensor on quadrature points
  ElementTypeMapArray<Real> conductivity_on_qpoints;

  /// vector \f[k \grad T\f] on quad points
  ElementTypeMapArray<Real> k_gradt_on_qpoints;

  /// external flux vector
  std::unique_ptr<Array<Real>> external_heat_rate;

  /// residuals array
  std::unique_ptr<Array<Real>> internal_heat_rate;

  /// boundary vector
  std::unique_ptr<Array<bool>> blocked_dofs;

  // realtime
  // Real time;

  /// capacity
  Real capacity;

  // conductivity matrix
  Matrix<Real> conductivity;

  // linear variation of the conductivity (for temperature dependent
  // conductivity)
  Real conductivity_variation;

  // reference temperature for the interpretation of temperature variation
  Real T_ref;

  // the biggest parameter of conductivity matrix
  // Real conductivitymax;

  bool need_to_reassemble_capacity{true};
  bool need_to_reassemble_capacity_lumped{true};
  Int temperature_release{0};
  Int conductivity_matrix_release{-1};
  std::unordered_map<GhostType, bool> initial_conductivity{{_not_ghost, true},
                                                           {_ghost, true}};
  std::unordered_map<GhostType, Int> conductivity_release{{_not_ghost, 0},
                                                          {_ghost, 0}};
};

} // namespace akantu

#endif /* AKANTU_HEAT_TRANSFER_MODEL_HH_ */
