/**
 * @file   heat_transfer_model.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Mon Mar 15 2021
 *
 * @brief  Model of Heat Transfer
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
#include "boundary_condition.hh"
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
                          public BoundaryCondition<HeatTransferModel>,
                          public DataAccessor<Element>,
                          public DataAccessor<UInt> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  HeatTransferModel(Mesh & mesh, UInt dim = _all_dimensions,
                    const ID & id = "heat_transfer_model",
                    ModelType model_type = ModelType::_heat_transfer_model,
                    std::shared_ptr<DOFManager> dof_manager = nullptr);

  ~HeatTransferModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// read one material file to instantiate all the materials
  virtual void readMaterials();

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
  virtual Real getStableTimeStep();

  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

public:
  /// calculate the lumped capacity vector for heat transfer problem
  void assembleCapacityLumped();

  /// compute the internal heat flux \todo Need code review
  virtual void assembleInternalHeatRate();

public:
  /// assemble the conductivity matrix
  virtual void assembleConductivityMatrix();

  /// assemble the conductivity matrix
  virtual void assembleCapacity();

  /// compute the capacity on quadrature points
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

  /// syncronize T or gradT directly by the model (for python wrapper)
  void synchronizeField(SynchronizationTag tag) { this->synchronize(tag); };

public:
  /// asign material properties to physical groups
  void assignPropertyToPhysicalGroup(const std::string & property_name,
                                     const std::string & group_name,
                                     Real value);
  /// asign material properties to physical groups
  void assignPropertyToPhysicalGroup(const std::string & property_name,
                                     const std::string & group_name,
                                     Matrix<Real> cond_matrix);

private:
  /// compute the conductivity tensor for each quadrature point in an array
  void computeConductivityOnQuadPoints(GhostType ghost_type);

  /// compute vector \f[k \grad T\f] for each quadrature point
  void computeKgradT(GhostType ghost_type);

  /// compute the thermal energy
  Real computeThermalEnergyByNode();

protected:
  /// calculate the lumped capacity vector for heat transfer problem (w
  /// ghost type)
  virtual void assembleCapacityLumped(GhostType ghost_type);
  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;
  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  inline UInt getNbData(const Array<UInt> & indexes,
                        const SynchronizationTag & tag) const override;
  inline void packData(CommunicationBuffer & buffer,
                       const Array<UInt> & indexes,
                       const SynchronizationTag & tag) const override;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<UInt> & indexes,
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
                       UInt spatial_dimension, ElementKind kind) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Density, density, Real);
  AKANTU_GET_MACRO(Capacity, capacity, Real);
  /// get the assembled heat flux
  AKANTU_GET_MACRO(InternalHeatRate, *internal_heat_rate, Array<Real> &);
  /// get the boundary vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool> &);
  /// get the external heat rate vector
  AKANTU_GET_MACRO(ExternalHeatRate, *external_heat_rate, Array<Real> &);
  /// get the temperature gradient
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureGradient,
                                         temperature_gradient, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConductivityOnQpoints,
                                         conductivity_on_qpoints, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureOnQpoints,
                                         temperature_on_qpoints, Real);
  /// get density array on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(DensityArray, density_array, Real);
  /// get capacity array on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(CapacityArray, capacity_array, Real);
  /// internal variables
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(KgradT, k_gradt_on_qpoints, Real);

  /// get the temperature
  AKANTU_GET_MACRO(Temperature, *temperature, Array<Real> &);
  /// get the temperature derivative
  AKANTU_GET_MACRO(TemperatureRate, *temperature_rate, Array<Real> &);
  /// get the temperature rate at qpoints
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureRateOnQpoints,
                                         temperature_rate_on_qpoints, Real);

  /// get the energy denominated by thermal
  Real getEnergy(const std::string & energy_id, ElementType type, UInt index);
  /// get the energy denominated by thermal
  Real getEnergy(const std::string & energy_id);

  /// get the thermal energy for a given element
  Real getThermalEnergy(ElementType type, UInt index);
  /// get the thermal energy for a given element
  Real getThermalEnergy();
  /// get the FEEngine object to apply Neumann BCs
  FEEngine & getFEEngineBoundary(const ID & name = "") override;
  /// evaluates need_to_reassemble_capacity boolean
  bool needToReassembleCapacity() const { return need_to_reassemble_capacity; }
  /// get the conductivity release
  UInt getConductivityRelease(GhostType ghost_type = _not_ghost) const {
    return conductivity_release.at(ghost_type);
  }
  /// get the conductivity matrix release
  AKANTU_GET_MACRO(ConductivityMatrixRelease, conductivity_matrix_release,
                   UInt);

protected:
  /* ----------------------------------------------------------------------- */
  template <class iterator>
  void getThermalEnergy(iterator Eth, Array<Real>::const_iterator<Real> T_it,
                        const Array<Real>::const_iterator<Real> & T_end,
                        Array<Real>::const_iterator<Real> capacity_it,
                        Array<Real>::const_iterator<Real> density_it) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// temperatures array
  std::unique_ptr<Array<Real>> temperature;

  /// temperatures derivatives array
  std::unique_ptr<Array<Real>> temperature_rate;

  /// increment array (@f$\delta \dot T@f$ or @f$\delta T@f$)
  std::unique_ptr<Array<Real>> increment;

  /// the density
  Real density;
  ElementTypeMapArray<Real> density_array;

  /// spatial temperature gradient
  ElementTypeMapArray<Real> temperature_gradient;

  /// temperature field on quadrature points
  ElementTypeMapArray<Real> temperature_on_qpoints;

  /// temperature rate field on quadrature points
  ElementTypeMapArray<Real> temperature_rate_on_qpoints;

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
  ElementTypeMapArray<Real> capacity_array;

  // conductivity matrix
  Matrix<Real> conductivity;
  ElementTypeMapArray<Real> initial_conductivity_array;

  // linear variation of the conductivity (for temperature dependent
  // conductivity)
  Real conductivity_variation;
  ElementTypeMapArray<Real> conductivity_variation_array;

  // reference temperature for the interpretation of temperature variation
  Real T_ref;
  ElementTypeMapArray<Real> T_ref_array;

  // the biggest parameter of conductivity matrix
  // Real conductivitymax;

  bool need_to_reassemble_capacity{true};
  bool need_to_reassemble_capacity_lumped{true};
  UInt temperature_release{0};
  UInt conductivity_matrix_release{UInt(-1)};
  std::unordered_map<GhostType, bool> initial_conductivity{{_not_ghost, true},
                                                           {_ghost, true}};
  std::unordered_map<GhostType, UInt> conductivity_release{{_not_ghost, 0},
                                                           {_ghost, 0}};
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "heat_transfer_model_inline_impl.hh"

#endif /* AKANTU_HEAT_TRANSFER_MODEL_HH_ */
