/**
 * @file   level_set_model.hh
 *
 * @author Danie Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Fri Dec 14 11:45:43 2012
 *
 * @brief  Level Set Method
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

#ifndef __AKANTU_LEVEL_SET_MODEL_HH__
#define __AKANTU_LEVEL_SET_MODEL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "model.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "dumpable.hh"
#include "solver.hh"
#include "geometry.hh"

namespace akantu {
  class IntegrationScheme1stOrder;
}

__BEGIN_AKANTU__

class LevelSetModel : public Model, public DataAccessor, public MeshEventHandler, public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange> MyFEEngineType;

  //LevelSetModel(UInt spatial_dimension,
  LevelSetModel(Mesh & mesh,
                UInt spatial_dimension = _all_dimensions,
                const ID & id = "level_set_model",
                const MemoryID & memory_id = 0);

  virtual ~LevelSetModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:

  /// generic function to initialize everything ready for explicit dynamics
  //void initFull(const std::string & material_file);
  void initFull();

  /// initialize the fem object of the boundary
  void initFEEngineBoundary(bool create_surface = false);

  /// allocate all vectors
  void initArrays();

  /// register the tags associated with the parallel synchronizer
  void initParallel(MeshPartition * partition, DataAccessor * data_accessor = NULL);

  /// initialize the model
  void initModel();

  /// init PBC synchronizer
  void initPBC();

  /// initialize the solver and the jacobian_matrix (called by initImplicit)
  void initSolver(SolverOptions & options = _solver_no_options);

  /// initialize the stuff for the implicit solver
  void initImplicit(SolverOptions & solver_options = _solver_no_options);

  void initPhi(geometry & geo, bool filter=false,  Real epsilon=0.0);

  /// function to print the contain of the class

  virtual void printself(__attribute__((unused)) std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {
  };

  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:

  /// compute and get the stable time step
  //Real getStableTimeStep();

  /// compute the rigth side had term: int(phi^{t-\Delta t}/\Delta t w )dV
  void updateRHS(bool reinit=false, Real Epsilon=0.0);

  /// Apply inflow boundary conditions
  //void ApplyInflowBC(bool matrix, bool reinit);

  void computeVReinit(Real Epsilon);

  /// calculate N_i N_j (w ghosttype)
  void assemblePhi(const GhostType & ghost_type, bool reinit=false);

  /// calculate v dot \grad(N_j) Ni (w ghosttype)
  //void assembleDiffusionPhi(const GhostType & ghost_type);

  /// calculate v dot \grad(N_j) Ni (w ghosttype)
  //void assemblePhi_t(const GhostType & ghost_type);

  bool testConvergenceResidual(Real tolerance, Real & norm);

  void solveStatic();

  /// calculate the lumped capacity vector for heat transfer problem
  //void assembleCapacityLumped(); @Daniel

  /// update the temperature from the temperature rate
  //void explicitPred();

  /// update the temperature rate from the increment
  //void explicitCorr();


  // /// initialize the heat flux
  // void initializeResidual(Array<Real> &temp);
  // /// initialize temperature
  // void initializeTemperature(Array<Real> &temp);

private:


  /// compute the heat flux on ghost types
  void updateRHS(const GhostType & ghost_type, bool reinit, Real Epsilon);

  void setIncrementFlagOn();

  void Shapes_SUPG(bool reinit=false);


  /// solve the system in temperature rate  @f$C\delta \dot T = q_{n+1} - C \dot T_{n}@f$
  //void solveExplicitLumped();

  /// calculate the lumped capacity vector for heat transfer problem (w ghosttype)
  //void assembleCapacityLumped(const GhostType & ghost_type);

  /// compute the conductivity tensor for each quadrature point in an array
  //void computeConductivityOnQuadPoints(const GhostType & ghost_type);

  /// compute vector k \grad T for each quadrature point
  //void computeKgradT(const GhostType & ghost_type);


  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:

  inline UInt getNbDataToPack(const Element & element,
                              SynchronizationTag tag) const;
  inline UInt getNbDataToUnpack(const Element & element,
                                SynchronizationTag tag) const;
  inline void packData(CommunicationBuffer & buffer,
                       const Element & element,
                       SynchronizationTag tag) const;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Element & element,
                         SynchronizationTag tag);

  inline UInt getNbDataToPack(SynchronizationTag tag) const;
  inline UInt getNbDataToUnpack(SynchronizationTag tag) const;
  inline void packData(CommunicationBuffer & buffer,
                       const UInt index,
                       SynchronizationTag tag) const;
  inline void unpackData(CommunicationBuffer & buffer,
                         const UInt index,
                         SynchronizationTag tag);

  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);

  void transportLevelSet(Array<Real> * velocity, Real delta_t=0.0, bool Assembly_A=false);
  void reinitializeLevelSet(Real delta_t, Real tol, UInt max_step, Real Epsilon);
  inline Real sign_phi(Real phi_q, Real epsilon=0.01){
    return phi_q/(sqrt(phi_q*phi_q+epsilon));
  };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  inline FEEngine & getFEEngineBoundary(std::string name = "");

  /// get the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// set the value of the time step
  AKANTU_SET_MACRO(TimeStep, time_step, Real);
  /// get the assembled heat flux
  AKANTU_GET_MACRO(Residual, *residual, Array<Real>&);
  /// get Stiffness matrix
  AKANTU_GET_MACRO(StiffnessMatrix, *stiffness_matrix, SparseMatrix &);
  /// get the lumped capacity
  //AKANTU_GET_MACRO(CapacityLumped, * capacity_lumped, Array<Real>&);
  /// get the boundary vector
  AKANTU_GET_MACRO(Boundary, * boundary, Array<bool>&);
  /// get the external flux vector
  //AKANTU_GET_MACRO(ExternalFlux, * external_flux, Array<Real>&);
  /// get the temperature gradient
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PhiGradient, phi_gradient, Real);
  /// get the conductivity on q points
  //AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConductivityOnQpoints, conductivity_on_qpoints, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PhiOnQpoints, phi_on_qpoints, Real);
  /// internal variables
  //AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(KGradtOnQpoints, k_gradt_on_qpoints, Real);
  //AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(IntBtKgT, int_bt_k_gT, Real);
  //AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(BtKgT, bt_k_gT, Real);
  /// get the temperature
  AKANTU_GET_MACRO(Phi, *phi, Array<Real> &);
  AKANTU_GET_MACRO(Velocity, *v, Array<Real> &);
  /// get the temperature derivative
  //AKANTU_GET_MACRO(TemperatureRate, *temperature_rate, Array<Real> &);
  /// get the equation number Array<Int>
  //AKANTU_GET_MACRO(EquationNumber, *equation_number, const Array<Int> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  Array<Real> * v;
protected:

  ///Level set function \phi
  Array<Real> * phi;

  /// increment of Phi
  Array<Real> * increment;

  /// flag defining if the increment must be computed or not
  bool increment_flag;

  ///Convection velocity


  ///\grad phi
  ElementTypeMapArray<Real> phi_gradient;

  ///\phi at quads points
  ElementTypeMapArray<Real> phi_on_qpoints;
  ElementTypeMapArray<Real> phi_on_qpoints_boundary;

  ///velocity at quads points
  ElementTypeMapArray<Real> v_on_qpoints;
  ElementTypeMapArray<Real> v_on_qpoints_boundary;
  ElementTypeMapArray<Real> v_r_on_qpoints;

  ///element by type
  ElementTypeMapArray<UInt> element_filter;
  ElementTypeMapArray<UInt> element_filter_boundary;

  ///SUPG
  ElementTypeMapArray<Real> shapes_SUPG;

  ///SUPG flag
  bool supg_flag;

  ///filtered flag
  bool filtered_flag;

  ///Filter parameter
  Real Filter_parameter;

  ///time step
  Real time_step;

  /// residuals array
  Array<Real> * residual;

  /// boundary vector
  Array<bool> * boundary;

  /// Matrix
  SparseMatrix * stiffness_matrix;

  /// jacobian matrix
  SparseMatrix * jacobian_matrix;

  /// Mesh
  Mesh & mesh;

  /// solver for implicit
  Solver * solver;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#include "level_set_model_inline_impl.cc"
#endif

/// standard output stream operator

inline std::ostream & operator <<(std::ostream & stream, const LevelSetModel & _this) {
  _this.printself(stream);
  return stream;
}


__END_AKANTU__



#endif /* __AKANTU_HEAT_TRANSFER_MODEL_HH__ */
