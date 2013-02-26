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

/* -------------------------------------------------------------------------- */
namespace akantu {
  class Material;
  class IntegrationScheme2ndOrder;
  class Contact;
  class SparseMatrix;
}

__BEGIN_AKANTU__

class SolidMechanicsModel : public Model, public DataAccessor, public MeshEventHandler, public Dumpable<DumperParaview> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  class NewMaterialElementsEvent : public NewElementsEvent {
  public:
    AKANTU_GET_MACRO_NOT_CONST(MaterialList, material, Vector<UInt> &);
    AKANTU_GET_MACRO(MaterialList, material, const Vector<UInt> &);
  protected:
    Vector<UInt> material;
  };

  typedef FEMTemplate<IntegratorGauss,ShapeLagrange> MyFEMType;

  SolidMechanicsModel(Mesh & mesh,
		      UInt spatial_dimension = 0,
		      const ID & id = "solid_mechanics_model",
		      const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize completely the model
  virtual void initFull(std::string material_file,
			AnalysisMethod method = _explicit_dynamic);

  /// initialize the fem object needed for boundary conditions
  void initFEMBoundary(bool create_surface = true);

  /// register the tags associated with the parallel synchronizer
  void initParallel(MeshPartition * partition, DataAccessor * data_accessor=NULL);

  /// allocate all vectors
  void initVectors();

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
  void initExplicit();

  /// initialize the array needed by updateResidual (residual, current_position)
  void initializeUpdateResidualData();

  /// update the current position vector
  void updateCurrentPosition();

  /// assemble the residual for the explicit scheme
  void updateResidual(bool need_initialize = true);

  /**
   * \brief compute the acceleration from the residual
   * this function is the explicit equivalent to solveDynamic in implicit
   * In the case of lumped mass just divide the residual by the mass
   * In the case of not lumped mass call solveDynamic<_acceleration_corrector>
   */
  void updateAcceleration();

  /// Solve the system @f[ A x = \alpha b @f] with A a lumped matrix
  void solveLumped(Vector<Real> & x, 
		   const Vector<Real> & A,
		   const Vector<Real> & b,
		   const Vector<bool> & boundary,
		   Real alpha);

  /// explicit integration predictor
  void explicitPred();

  /// explicit integration corrector
  void explicitCorr();

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

  /// solve @f[ A\delta u = f_ext - f_int @f] in displacement
  void solveDynamic();

  /// solve Ku = f
  void solveStatic();

  /// test the convergence (norm of increment)
  bool testConvergenceIncrement(Real tolerance);
  bool testConvergenceIncrement(Real tolerance, Real & error);

  /// test the convergence (norm of residual)
  bool testConvergenceResidual(Real tolerance);
  bool testConvergenceResidual(Real tolerance, Real & error);

  /// implicit time integration predictor
  void implicitPred();

  /// implicit time integration corrector
  void implicitCorr();

protected:
  /// finish the computation of residual to solve in increment
  void updateResidualInternal();

  /// compute the support reaction and store it in force
  void updateSupportReaction();

  /// compute A and solve @f[ A\delta u = f_ext - f_int @f]
  template<NewmarkBeta::IntegrationSchemeCorrectorType type>
  void solveDynamic(Vector<Real> & increment);

  /* ------------------------------------------------------------------------ */
  /* Explicit/Implicit                                                        */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the stresses
  void computeStresses();

  /* ------------------------------------------------------------------------ */
  /* Boundaries (solid_mechanics_model_boundary.cc)                           */
  /* ------------------------------------------------------------------------ */
public:
  class SurfaceLoadFunctor {
  public:
    virtual void traction(__attribute__ ((unused)) const types::Vector<Real> & position,
			  __attribute__ ((unused)) types::Vector<Real> & traction,
			  __attribute__ ((unused)) const types::Vector<Real> & normal,
			  __attribute__ ((unused)) Surface surface_id) {
      AKANTU_DEBUG_TO_IMPLEMENT();
    }

    virtual void stress(__attribute__ ((unused)) const types::Vector<Real> & position,
			__attribute__ ((unused)) types::RMatrix & stress,
			__attribute__ ((unused)) const types::Vector<Real> & normal,
			__attribute__ ((unused)) Surface surface_id) {
      AKANTU_DEBUG_TO_IMPLEMENT();
    }
  };

  /// compute force vector from a function(x,y,z) that describe stresses
  void computeForcesFromFunction(BoundaryFunction in_function,
				 BoundaryFunctionType function_type) __attribute__((deprecated));

  template<class Functor>
  void computeForcesFromFunction(Functor & functor, BoundaryFunctionType function_type);

  /// integrate a force on the boundary by providing a stress tensor
  void computeForcesByStressTensor(const Vector<Real> & stresses,
				   const ElementType & type,
				   const GhostType & ghost_type);

  /// integrate a force on the boundary by providing a traction vector
  void computeForcesByTractionVector(const Vector<Real> & tractions,
				     const ElementType & type,
				     const GhostType & ghost_type);

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

  /// assemble the mass matrix
  void assembleMass();


protected:
  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumped(GhostType ghost_type);

  void assembleMass(GhostType ghost_type);

  /// fill a vector of rho
  void computeRho(Vector<Real> & rho,
		  ElementType type,
		  GhostType ghost_type);


  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:

  inline virtual UInt getNbDataForElements(const Vector<Element> & elements,
					   SynchronizationTag tag) const;

  inline virtual void packElementData(CommunicationBuffer & buffer,
				      const Vector<Element> & elements,
				      SynchronizationTag tag) const;

  inline virtual void unpackElementData(CommunicationBuffer & buffer,
					const Vector<Element> & elements,
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
  inline void splitElementByMaterial(const Vector<Element> & elements,
				     Vector<Element> * elements_per_mat) const;

  inline virtual void packBarycenter(CommunicationBuffer & buffer,
				     const Vector<Element> & elements) const;

  inline virtual void unpackBarycenter(CommunicationBuffer & buffer,
				       const Vector<Element> & elements,
				       SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  virtual void onNodesAdded  (const Vector<UInt> & nodes_list,
			      const NewNodesEvent & event);
  virtual void onNodesRemoved(const Vector<UInt> & element_list,
			      const Vector<UInt> & new_numbering,
			      const RemovedNodesEvent & event);
  virtual void onElementsAdded  (const Vector<Element> & nodes_list,
				 const NewElementsEvent & event);
  virtual void onElementsRemoved(const Vector<Element> & element_list,
				 const ByElementTypeUInt & new_numbering,
				 const RemovedElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpField(const std::string & field_id);
  virtual void addDumpFieldVector(const std::string & field_id);
  virtual void addDumpFieldTensor(const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// set the value of the time step
  AKANTU_SET_MACRO(TimeStep, time_step, Real);

  /// get the value of the conversion from forces/ mass to acceleration
  AKANTU_GET_MACRO(F_M2A, f_m2a, Real);
  /// set the value of the conversion from forces/ mass to acceleration
  AKANTU_SET_MACRO(F_M2A, f_m2a, Real);

  /// get the SolidMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement,    *displacement,           Vector<Real> &);
  /// get the SolidMechanicsModel::current_position vector \warn only consistent
  /// after a call to SolidMechanicsModel::updateCurrentPosition
  AKANTU_GET_MACRO(CurrentPosition, *current_position, const Vector<Real> &);
  /// get  the SolidMechanicsModel::increment  vector \warn  only  consistent if
  /// SolidMechanicsModel::setIncrementFlagOn has been called before
  AKANTU_GET_MACRO(Increment,       *increment,              Vector<Real> &);
  /// get the lumped SolidMechanicsModel::mass vector
  AKANTU_GET_MACRO(Mass,            *mass,                   Vector<Real> &);
  /// get the SolidMechanicsModel::velocity vector
  AKANTU_GET_MACRO(Velocity,        *velocity,               Vector<Real> &);
  /// get    the    SolidMechanicsModel::acceleration    vector,   updated    by
  /// SolidMechanicsModel::updateAcceleration
  AKANTU_GET_MACRO(Acceleration,    *acceleration,           Vector<Real> &);
  /// get the SolidMechanicsModel::force vector (boundary forces)
  AKANTU_GET_MACRO(Force,           *force,                  Vector<Real> &);
  /// get     the    SolidMechanicsModel::residual    vector,     computed    by
  /// SolidMechanicsModel::updateResidual
  AKANTU_GET_MACRO(Residual,        *residual,               Vector<Real> &);
  /// get the SolidMechanicsModel::boundary vector
  AKANTU_GET_MACRO(Boundary,        *boundary,               Vector<bool> &);
  /// get the SolidMechnicsModel::incrementAcceleration vector
  AKANTU_GET_MACRO(IncrementAcceleration, *increment_acceleration, Vector<Real> &);

  /// get the value of the SolidMechanicsModel::increment_flag
  AKANTU_GET_MACRO(IncrementFlag, increment_flag, bool);

  /// get  the SolidMechanicsModel::element_material  vector corresponding  to a
  /// given akantu::ElementType
  //  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementMaterial, element_material, UInt);
  //  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt);
  //AKANTU_GET_MACRO(ElementMaterial, element_material, const ByElementTypeVector<UInt> &);

  /// vectors containing local material element index for each global element index
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementIndexByMaterial, element_index_by_material, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementIndexByMaterial, element_index_by_material, UInt);
  AKANTU_GET_MACRO(ElementIndexByMaterial, element_index_by_material, const ByElementTypeVector<UInt> &);

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

  /// compute the external work (for impose displacement, the velocity should be given too)
  Real getExternalWork();

  /// get the energies
  Real getEnergy(const std::string & energy_id);

  /// compute the energy for energy
  Real getEnergy(const std::string & energy_id, ElementType & type, UInt index);

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

  inline FEM & getFEMBoundary(const ID & name = "");

  /// get integrator
  AKANTU_GET_MACRO(Integrator, *integrator, const IntegrationScheme2ndOrder &);

  /// get access to the internal solver
  AKANTU_GET_MACRO(Solver, *solver, Solver &);

protected:
  /// compute the stable time step
  Real getStableTimeStep(const GhostType & ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// time step
  Real time_step;

  /// conversion coefficient form force/mass to acceleration
  Real f_m2a;

  /// displacements array
  Vector<Real> * displacement;

  /// lumped mass array
  Vector<Real> * mass;

  /// velocities array
  Vector<Real> * velocity;

  /// accelerations array
  Vector<Real> * acceleration;

  /// accelerations array
  Vector<Real> * increment_acceleration;

  /// forces array
  Vector<Real> * force;

  /// residuals array
  Vector<Real> * residual;

  /// boundaries array
  Vector<bool> * boundary;

  /// array of current position used during update residual
  Vector<Real> * current_position;

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
  Vector<Real> * increment;

  /// flag defining if the increment must be computed or not
  bool increment_flag;

  /// solver for implicit
  Solver * solver;

  /// object to resolve the contact
  Contact * contact;

  /// the spatial dimension
  UInt spatial_dimension;

  AnalysisMethod method;

  Synchronizer * synch_parallel;
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
