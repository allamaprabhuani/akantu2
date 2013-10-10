/**
 * @file   solid_mechanics_model.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Jul 27 18:15:37 2010
 *
 * @brief  Model of Solid Mechanics
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

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_HH__


/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "model.hh"
#include "data_accessor.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "aka_types.hh"
#include "integration_scheme_2nd_order.hh"
#include "solver.hh"
#include "dumpable.hh"
#include "boundary_condition.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
  class Material;
  class IntegrationScheme2ndOrder;
  class Contact;
  class SparseMatrix;
}

__BEGIN_AKANTU__

class SolidMechanicsModel : public Model, public DataAccessor,
                            public MeshEventHandler,
                            public Dumpable,
                            public BoundaryCondition<SolidMechanicsModel> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  class NewMaterialElementsEvent : public NewElementsEvent {
  public:
    AKANTU_GET_MACRO_NOT_CONST(MaterialList, material, Array<UInt> &);
    AKANTU_GET_MACRO(MaterialList, material, const Array<UInt> &);
  protected:
    Array<UInt> material;
  };

  typedef FEMTemplate<IntegratorGauss,ShapeLagrange> MyFEMType;

  SolidMechanicsModel(Mesh & mesh,
                      UInt spatial_dimension = _all_dimensions,
                      const ID & id = "solid_mechanics_model",
                      const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize completely the model
  virtual void initFull(std::string material_file,
                        AnalysisMethod method = _explicit_lumped_mass);

  /// initialize the fem object needed for boundary conditions
  void initFEMBoundary();

  /// register the tags associated with the parallel synchronizer
  void initParallel(MeshPartition * partition, DataAccessor * data_accessor=NULL);

  /// allocate all vectors
  void initArrays();

  /// allocate all vectors
  void initArraysPreviousDisplacment();

  /// initialize all internal arrays for materials
  void initMaterials();

  /// initialize the model
  virtual void initModel();

  /// init PBC synchronizer
  void initPBC();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* PBC                                                                      */
  /* ------------------------------------------------------------------------ */
public:

  /// change the equation number for proper assembly when using PBC
  void changeEquationNumberforPBC(std::map<UInt,UInt> & pbc_pair);

  /// synchronize Residual for output
  void synchronizeResidual();

protected:

  /// register PBC synchronizer
  void registerPBCSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Explicit                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the stuff for the explicit scheme
  void initExplicit(AnalysisMethod analysis_method = _explicit_lumped_mass);

  bool isExplicit()
  { return method == _explicit_lumped_mass || method == _explicit_consistent_mass; }

  /// initialize the array needed by updateResidual (residual, current_position)
  void initializeUpdateResidualData();

  /// update the current position vector
  void updateCurrentPosition();

  /// assemble the residual for the explicit scheme
  virtual void updateResidual(bool need_initialize = true);

  /**
   * \brief compute the acceleration from the residual
   * this function is the explicit equivalent to solveDynamic in implicit
   * In the case of lumped mass just divide the residual by the mass
   * In the case of not lumped mass call solveDynamic<_acceleration_corrector>
   */
  void updateAcceleration();

  void updateIncrement();
  void updatePreviousDisplacement();
 
  /// Solve the system @f[ A x = \alpha b @f] with A a lumped matrix
  void solveLumped(Array<Real> & x,
                   const Array<Real> & A,
                   const Array<Real> & b,
                   const Array<bool> & boundary,
                   Real alpha);

  /// explicit integration predictor
  void explicitPred();

  /// explicit integration corrector
  void explicitCorr();

public:
  void solveStep();

  /* ------------------------------------------------------------------------ */
  /* Implicit                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the solver and the jacobian_matrix (called by initImplicit)
  void initSolver(SolverOptions & options = _solver_no_options);

  /// initialize the stuff for the implicit solver
  void initImplicit(bool dynamic = false,
                    SolverOptions & solver_options = _solver_no_options);

  /// solve Ma = f to get the initial acceleration
  void initialAcceleration();

  /// assemble the stiffness matrix
  void assembleStiffnessMatrix();

public:
  template<SolveConvergenceMethod method, SolveConvergenceCriteria criteria>
  void solveStep(Real tolerance,
                 UInt max_iteration = 100);

public:
  /// solve @f[ A\delta u = f_ext - f_int @f] in displacement
  void solveDynamic() __attribute__((deprecated));

  /// solve Ku = f
  void solveStatic() __attribute__((deprecated));

  /// solve Ku = f
  void solveStatic(Array<bool> & boundary_normal, Array<Real> & EulerAngles);

  /// test if the system is converged
  template<SolveConvergenceCriteria criteria>
  bool testConvergence(Real tolerance, Real & error);

  /// test the convergence (norm of increment)
  bool testConvergenceIncrement(Real tolerance) __attribute__((deprecated));
  bool testConvergenceIncrement(Real tolerance, Real & error) __attribute__((deprecated));

  /// test the convergence (norm of residual)
  bool testConvergenceResidual(Real tolerance) __attribute__((deprecated));
  bool testConvergenceResidual(Real tolerance, Real & error) __attribute__((deprecated));

  /// create and return the velocity damping matrix
  SparseMatrix & initVelocityDampingMatrix();

  /// implicit time integration predictor
  void implicitPred();

  /// implicit time integration corrector
  void implicitCorr();

protected:
  /// finish the computation of residual to solve in increment
  void updateResidualInternal();

  /// compute the support reaction and store it in force
  void updateSupportReaction();

public://protected: Daniel changed it just for a test
  /// compute A and solve @f[ A\delta u = f_ext - f_int @f]
  template<NewmarkBeta::IntegrationSchemeCorrectorType type>
  void solve(Array<Real> & increment);

  /* ------------------------------------------------------------------------ */
  /* Explicit/Implicit                                                        */
  /* ------------------------------------------------------------------------ */
public:
  /// Update the stresses for the computation of the residual of the Stiffness
  /// matrix in the case of finite deformation
  void computeStresses();

  /// synchronize the ghost element boundaries values
  void synchronizeBoundaries();

  /* ------------------------------------------------------------------------ */
  /* Materials (solid_mechanics_model_material.cc)                            */
  /* ------------------------------------------------------------------------ */
public:
  /// register a material in the dynamic database
  Material & registerNewMaterial(const ID & mat_type,
                                 const std::string & opt_param = "");

  template <typename M>
  Material & registerNewCustomMaterial(const ID & mat_type,
                                       const std::string & opt_param = "");
  /// read the material files to instantiate all the materials
  void readMaterials(const std::string & filename);

  /// read a custom material with a keyword and class as template
  template <typename M>
  UInt readCustomMaterial(const std::string & filename,
                                const std::string & keyword);

  /// Use a UIntData in the mesh to specify the material to use per element
  void setMaterialIDsFromIntData(const std::string & data_name);

  /* ------------------------------------------------------------------------ */
  /* Mass (solid_mechanics_model_mass.cc)                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// assemble the lumped mass matrix
  void assembleMassLumped();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMass();


protected:
  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumped(GhostType ghost_type);

  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMass(GhostType ghost_type);

  /// fill a vector of rho
  void computeRho(Array<Real> & rho,
                  ElementType type,
                  GhostType ghost_type);


  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:

  inline virtual UInt getNbDataForElements(const Array<Element> & elements,
                                           SynchronizationTag tag) const;

  inline virtual void packElementData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      SynchronizationTag tag) const;

  inline virtual void unpackElementData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        SynchronizationTag tag);

  inline virtual UInt getNbDataToPack(SynchronizationTag tag) const;
  inline virtual UInt getNbDataToUnpack(SynchronizationTag tag) const;

  inline virtual void packData(CommunicationBuffer & buffer,
                               const UInt index,
                               SynchronizationTag tag) const;

  inline virtual void unpackData(CommunicationBuffer & buffer,
                                 const UInt index,
                                 SynchronizationTag tag);

protected:
  inline void splitElementByMaterial(const Array<Element> & elements,
                                     Array<Element> * elements_per_mat) const;

  inline virtual void packBarycenter(const Mesh & mesh,
                                     CommunicationBuffer & buffer,
                                     const Array<Element> & elements,
                                     SynchronizationTag tag) const;

  inline virtual void unpackBarycenter(const Mesh & mesh,
                                       CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  virtual void onNodesAdded  (const Array<UInt> & nodes_list,
                              const NewNodesEvent & event);
  virtual void onNodesRemoved(const Array<UInt> & element_list,
                              const Array<UInt> & new_numbering,
                              const RemovedNodesEvent & event);
  virtual void onElementsAdded  (const Array<Element> & nodes_list,
                                 const NewElementsEvent & event);
  virtual void onElementsRemoved(const Array<Element> & element_list,
                                 const ByElementTypeUInt & new_numbering,
                                 const RemovedElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);

  virtual void addDumpBoundaryField(const std::string & field_id,
				    const std::string & boundary_name);
  virtual void addDumpBoundaryFieldToDumper(const std::string & dumper_name,
					    const std::string & field_id,
					    const std::string & boundary_name);
  virtual void removeDumpBoundaryField(const std::string & field_id,
				       const std::string & boundary_name);
  virtual void removeDumpBoundaryFieldFromDumper(const std::string & dumper_name,
						 const std::string & field_id,
						 const std::string & boundary_name);

  virtual void addDumpFieldVectorToDumper(const std::string & dumper_name,
					  const std::string & field_id);

  virtual void addDumpBoundaryFieldVector(const std::string & field_id,
					  const std::string & boundary_name);
  virtual void addDumpBoundaryFieldVectorToDumper(const std::string & dumper_name,
						  const std::string & field_id,
						  const std::string & boundary_name);

  virtual void addDumpFieldTensorToDumper(const std::string & dumper_name,
					  const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// set the value of the time step
  void setTimeStep(Real time_step);

  /// get the value of the conversion from forces/ mass to acceleration
  AKANTU_GET_MACRO(F_M2A, f_m2a, Real);
  /// set the value of the conversion from forces/ mass to acceleration
  AKANTU_SET_MACRO(F_M2A, f_m2a, Real);

  /// get the SolidMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement,    *displacement,           Array<Real> &);
  /// get the SolidMechanicsModel::previous_displacement vector
  AKANTU_GET_MACRO(PreviousDisplacement, *previous_displacement, Array<Real> &);
  /// get the SolidMechanicsModel::current_position vector \warn only consistent
  /// after a call to SolidMechanicsModel::updateCurrentPosition
  AKANTU_GET_MACRO(CurrentPosition, *current_position, const Array<Real> &);
  /// get  the SolidMechanicsModel::increment  vector \warn  only  consistent if
  /// SolidMechanicsModel::setIncrementFlagOn has been called before
  AKANTU_GET_MACRO(Increment,       *increment,              Array<Real> &);
  /// get the lumped SolidMechanicsModel::mass vector
  AKANTU_GET_MACRO(Mass,            *mass,                   Array<Real> &);
  /// get the SolidMechanicsModel::velocity vector
  AKANTU_GET_MACRO(Velocity,        *velocity,               Array<Real> &);
  /// get    the    SolidMechanicsModel::acceleration    vector,   updated    by
  /// SolidMechanicsModel::updateAcceleration
  AKANTU_GET_MACRO(Acceleration,    *acceleration,           Array<Real> &);
  /// get the SolidMechanicsModel::force vector (boundary forces)
  AKANTU_GET_MACRO(Force,           *force,                  Array<Real> &);
  /// get     the    SolidMechanicsModel::residual    vector,     computed    by
  /// SolidMechanicsModel::updateResidual
  AKANTU_GET_MACRO(Residual,        *residual,               Array<Real> &);
  /// get the SolidMechanicsModel::boundary vector
  AKANTU_GET_MACRO(Boundary,        *boundary,               Array<bool> &);
  /// get the SolidMechnicsModel::incrementAcceleration vector
  AKANTU_GET_MACRO(IncrementAcceleration, *increment_acceleration, Array<Real> &);

  /// get the value of the SolidMechanicsModel::increment_flag
  AKANTU_GET_MACRO(IncrementFlag, increment_flag, bool);

  /// get  the SolidMechanicsModel::element_material  vector corresponding  to a
  /// given akantu::ElementType
  //  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementMaterial, element_material, UInt);
  //  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt);
  //AKANTU_GET_MACRO(ElementMaterial, element_material, const ByElementTypeArray<UInt> &);

  /// vectors containing local material element index for each global element index
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementIndexByMaterial, element_index_by_material, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementIndexByMaterial, element_index_by_material, UInt);
  AKANTU_GET_MACRO(ElementIndexByMaterial, element_index_by_material, const ByElementTypeArray<UInt> &);

  /// get a particular material
  inline Material & getMaterial(UInt mat_index);
  inline const Material & getMaterial(UInt mat_index) const;

  /// give the number of materials
  inline UInt getNbMaterials() const { return materials.size(); };

  /// give the material internal index from its id
  Int getInternalIndexFromID(const ID & id) const;

  /// compute the stable time step
  Real getStableTimeStep();

  /// compute the potential energy
  Real getPotentialEnergy();

  /// compute the kinetic energy
  Real getKineticEnergy();
  Real getKineticEnergy(const ElementType & type, UInt index);

  /// compute the external work (for impose displacement, the velocity should be given too)
  Real getExternalWork();

  /// get the energies
  Real getEnergy(const std::string & energy_id);

  /// compute the energy for energy
  Real getEnergy(const std::string & energy_id, const ElementType & type, UInt index);

  /// set the Contact object
  AKANTU_SET_MACRO(Contact, contact, Contact *);

  /**
   * @brief set the SolidMechanicsModel::increment_flag  to on, the activate the
   * update of the SolidMechanicsModel::increment vector
   */
  void setIncrementFlagOn();

  /// get the stiffness matrix
  AKANTU_GET_MACRO(StiffnessMatrix, *stiffness_matrix, SparseMatrix &);

  /// get the mass matrix
  AKANTU_GET_MACRO(MassMatrix, *mass_matrix, SparseMatrix &);

  /// get the velocity damping matrix
  AKANTU_GET_MACRO(VelocityDampingMatrix, *velocity_damping_matrix, SparseMatrix &);

  /// get the FEM object to integrate or interpolate on the boundary
  inline FEM & getFEMBoundary(const ID & name = "");

  /// get integrator
  AKANTU_GET_MACRO(Integrator, *integrator, const IntegrationScheme2ndOrder &);

  /// get access to the internal solver
  AKANTU_GET_MACRO(Solver, *solver, Solver &);

  /// get synchronizer
  AKANTU_GET_MACRO(Synchronizer, *synch_parallel, const DistributedSynchronizer &);

protected:
  /// compute the stable time step
  Real getStableTimeStep(const GhostType & ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  friend class BoundaryCondition;

  /// time step
  Real time_step;

  /// conversion coefficient form force/mass to acceleration
  Real f_m2a;

  /// displacements array
  Array<Real> * displacement;

  /// displacements array at the previous time step (used in finite deformation)
  Array<Real> * previous_displacement;

  /// lumped mass array
  Array<Real> * mass;

  /// velocities array
  Array<Real> * velocity;

  /// accelerations array
  Array<Real> * acceleration;

  /// accelerations array
  Array<Real> * increment_acceleration;

  /// forces array
  Array<Real> * force;

  /// residuals array
  Array<Real> * residual;

  /// boundaries array
  Array<bool> * boundary;

  /// array of current position used during update residual
  Array<Real> * current_position;

  /// mass matrix
  SparseMatrix * mass_matrix;

  /// velocity damping matrix
  SparseMatrix * velocity_damping_matrix;

  /// stiffness matrix
  SparseMatrix * stiffness_matrix;

  /// jacobian matrix @f[A = cM + dD + K@f] with @f[c = \frac{1}{\beta \Delta t^2}, d = \frac{\gamma}{\beta \Delta t} @f]
  SparseMatrix * jacobian_matrix;

  /// vectors containing local material element index for each global element index
  ByElementTypeUInt element_index_by_material;

  /// list of used materials
  std::vector<Material *> materials;

  /// integration scheme of second order used
  IntegrationScheme2ndOrder * integrator;

  /// increment of displacement
  Array<Real> * increment;

  /// flag defining if the increment must be computed or not
  bool increment_flag;

  /// solver for implicit
  Solver * solver;

  /// object to resolve the contact
  Contact * contact;

  AnalysisMethod method;

  DistributedSynchronizer * synch_parallel;
};

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "parser.hh"
#include "material.hh"

__BEGIN_AKANTU__

#include "solid_mechanics_model_tmpl.hh"

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "solid_mechanics_model_inline_impl.cc"
#endif

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SolidMechanicsModel & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_HH__ */
