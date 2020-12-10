/**
 * @file   fluid_diffusion_model.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Aug 2019
 *
 * @brief  Transient fluid diffusion model based on the Heat Transfer Model
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
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "fluid_diffusion_model_event_handler.hh"
#include "model.hh"
#if defined(AKANTU_COHESIVE_ELEMENT)
#include "material_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
#endif

/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_FLUID_DIFFUSION_MODEL_HH__
#define __AKANTU_FLUID_DIFFUSION_MODEL_HH__

namespace akantu {
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

namespace akantu {

class FluidDiffusionModel
    : public Model,
      public DataAccessor<Element>,
      public DataAccessor<UInt>,
      public EventHandlerManager<FluidDiffusionModelEventHandler> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  FluidDiffusionModel(Mesh & mesh, UInt dim = _all_dimensions,
                      const ID & id = "fluid_diffusion_model",
                      const MemoryID & memory_id = 0);

    ~FluidDiffusionModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// read one material file to instantiate all the materials
  void readMaterials();

  /// allocate all vectors
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /// initialize the model
  void initModel() override;

  void predictor() override;

  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool is_converged = false) override;

  /// compute the heat flux
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID &) override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID &) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID &) override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  /// resize all fields on nodes added event
  void resizeFields();

public:
  /// get coordinates of qpoints in same order as other elemental fields
  Array<Real> getQpointsCoord();

#ifdef AKANTU_COHESIVE_ELEMENT
  /// get apperture at nodes from displacement field of cohesive model
  void getApertureOnQpointsFromCohesive(
      const SolidMechanicsModelCohesive & coh_model, bool first_time = false);

  /// insert fluid elements based on cohesive elements and their damage
  void updateFluidElementsFromCohesive(
      const SolidMechanicsModelCohesive & coh_model);
#endif

  void applyExternalFluxAtElementGroup(const Real & rate,
                                       const ElementGroup & source_facets,
                                       GhostType ghost_type = _not_ghost);

  void injectIntoFacetsByCoord(const Vector<Real> & position,
                               const Real & injection_rate);
  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// compute and get the stable time step
  Real getStableTimeStep();

  /// set the stable timestep
  void setTimeStep(Real dt, const ID & solver_id = "") override;

  // temporary protection to prevent bad usage: should check for bug
protected:
  /// compute the internal heat flux \todo Need code review: currently not
  /// public method
  void assembleInternalFlux();

public:
  /// calculate the lumped capacity vector for heat transfer problem
  void assembleCapacityLumped();

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  void onElementsAdded(const Array<Element> & element_list,
                       const NewElementsEvent & /*event*/) override;
  /* ------------------------------------------------------------------------ */
  /* Methods for static                                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// assemble the conductivity matrix
  void assemblePermeabilityMatrix();

  /// assemble the conductivity matrix
  void assembleCapacity();

  /// compute the capacity on quadrature points
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

private:
  /// calculate the lumped capacity vector for heat transfer problem (w
  /// ghost type)
  void assembleCapacityLumped(const GhostType & ghost_type);

  /// assemble the conductivity matrix (w/ ghost type)
  void assemblePermeabilityMatrix(const GhostType & ghost_type);

  /// compute the conductivity tensor for each quadrature point in an array
  void computePermeabilityOnQuadPoints(const GhostType & ghost_type);

  /// compute vector k \grad T for each quadrature point
  void computeKgradP(const GhostType & ghost_type);

  /// compute the thermal energy
  Real computeThermalEnergyByNode();

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

  virtual void dump(const std::string & dumper_name);
  virtual void dump(const std::string & dumper_name, UInt step);
  virtual void dump(const std::string & dumper_name, Real time, UInt step);

  void dump() override;

  virtual void dump(UInt step);

  virtual void dump(Real time, UInt step);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Compressibility, compressibility, Real);
  AKANTU_GET_MACRO(Pushability, pushability, Real);
  /// get the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// get the assembled heat flux
  AKANTU_GET_MACRO(InternalFlux, *internal_flux, Array<Real> &);
  /// get the boundary vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool> &);
  /// get the external heat rate vector
  AKANTU_GET_MACRO(ExternalFlux, *external_flux, Array<Real> &);
  /// get the temperature gradient
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PressureGradient, pressure_gradient,
                                         Real);
  /// get the permeability on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PermeabilityOnQpoints,
                                         permeability_on_qpoints, Real);
  /// get the pressure on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PressureOnQpoints, pressure_on_qpoints,
                                         Real);
  /// get the pressure change on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(DeltaPressureOnQpoints,
                                         delta_pres_on_qpoints, Real);
  /// get the aperture on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ApertureOnQpoints, aperture_on_qpoints,
                                         Real);
  /// get the aperture on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ApertureOnQpoints, aperture_on_qpoints,
                                   Real);
  /// internal variables
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(KgradP, k_gradp_on_qpoints, Real);
  /// get the temperature
  AKANTU_GET_MACRO(Pressure, *pressure, Array<Real> &);
  /// get the temperature derivative
  AKANTU_GET_MACRO(PressureRate, *pressure_rate, Array<Real> &);

  inline bool isModelVelocityDependent() const { return use_aperture_speed; }
  // /// get the energy denominated by thermal
  // Real getEnergy(const std::string & energy_id, const ElementType & type,
  //                UInt index);
  // /// get the energy denominated by thermal
  // Real getEnergy(const std::string & energy_id);

  // /// get the thermal energy for a given element
  // Real getThermalEnergy(const ElementType & type, UInt index);
  // /// get the thermal energy for a given element
  // Real getThermalEnergy();

protected:
  /* ------------------------------------------------------------------------ */
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ----------------------------------------------------------------------- */
  // template <class iterator>
  // void getThermalEnergy(iterator Eth, Array<Real>::const_iterator<Real> T_it,
  //                       Array<Real>::const_iterator<Real> T_end) const;

  template <typename T>
  void allocNodalField(Array<T> *& array, const ID & name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// number of iterations
  UInt n_iter;

  /// time step
  Real time_step;

  /// pressure array
  Array<Real> * pressure{nullptr};

  /// pressure derivatives array
  Array<Real> * pressure_rate{nullptr};

  // /// increment array (@f$\delta \dot P@f$ or @f$\delta P@f$)
  // Array<Real> * increment{nullptr};

  /// the speed of the changing temperature
  ElementTypeMapArray<Real> pressure_gradient;

  /// pressure field on quadrature points
  ElementTypeMapArray<Real> pressure_on_qpoints;

  /// pressure change on quadrature points
  ElementTypeMapArray<Real> delta_pres_on_qpoints;

  /// conductivity tensor on quadrature points
  ElementTypeMapArray<Real> permeability_on_qpoints;

  /// aperture on quadrature points
  ElementTypeMapArray<Real> aperture_on_qpoints;

  /// aperture on quadrature points
  ElementTypeMapArray<Real> prev_aperture_on_qpoints;

  /// vector k \grad T on quad points
  ElementTypeMapArray<Real> k_gradp_on_qpoints;

  /// external flux vector
  Array<Real> * external_flux{nullptr};

  /// residuals array
  Array<Real> * internal_flux{nullptr};

  /// boundary vector
  Array<bool> * blocked_dofs{nullptr};

  // viscosity
  Real viscosity{0.};

  /// compressibility Cf = 1/pho drho/dP
  Real compressibility{0.};

  // default value of aperture on a newly added nodesynchronizer
  Real default_aperture{0.};

  bool need_to_reassemble_capacity{true};
  bool need_to_reassemble_capacity_lumped{true};
  UInt pressure_release{0};
  UInt solution_release{0};
  UInt permeability_matrix_release{0};

  std::unordered_map<GhostType, bool> initial_permeability{{_not_ghost, true},
                                                           {_ghost, true}};
  std::unordered_map<GhostType, UInt> permeability_release{{_not_ghost, 0},
                                                           {_ghost, 0}};
  bool use_aperture_speed{false};

  // pushability Ch = 1/h dh/dP
  Real pushability{0.};

  // damage at which cohesive element will be duplicated in fluid mesh
  Real insertion_damage{0.};
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "fluid_diffusion_model_inline_impl.hh"

#endif /* __AKANTU_FLUID_DIFFUSION_MODEL_HH__ */
