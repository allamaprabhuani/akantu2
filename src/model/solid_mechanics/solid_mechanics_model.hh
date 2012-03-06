/**
 * @file   solid_mechanics_model.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date[creation]            Thu Jul 22 11:51:06 2010
 * @date[last modification]   Thu Oct 14 14:00:06 2010
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
#include "material.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "integrator_cohesive.hh"
#include "shape_cohesive.hh"
#include "aka_types.hh"
#include "integration_scheme_2nd_order.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  //  class Material;
  class IntegrationScheme2ndOrder;
  class Contact;
  class Solver;
  class SparseMatrix;
}

__BEGIN_AKANTU__

class SolidMechanicsModel : public Model, public DataAccessor {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEMTemplate<IntegratorGauss,ShapeLagrange> MyFEMType;
  typedef FEMTemplate< IntegratorCohesive<IntegratorGauss>, ShapeCohesive<ShapeLagrange> > MyFEMCohesiveType;

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
  void initFull(std::string material_file = "",
                bool implicit_scheme = false,
                bool implicit_dynamic = false);

  /// initialize the fem object needed for boundary conditions
  void initFEMBoundary(bool create_surface = true);

  /// register the tags associated with the parallel synchronizer
  void initParallel(MeshPartition * partition, DataAccessor * data_accessor=NULL);

  /// allocate all vectors
  void initVectors();

  /// initialize all internal arrays for materials
  void initMaterials();

  /// initialize the model
  void initModel();

  /// init PBC synchronizer
  void initPBC();

  /// initialize the structures for cohesive elements
  void initCohesive();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* PBC                                                                      */
  /* ------------------------------------------------------------------------ */
public:

  /// change the equation number for proper assembly when using PBC
  void changeEquationNumberforPBC(std::map<UInt,UInt> & pbc_pair);

private:

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

  /// explicit integration predictor
  void explicitPred();

  /// explicit integration corrector
  void explicitCorr();

  /* ------------------------------------------------------------------------ */
  /* Implicit                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the solver and the jacobian_matrix (called by initImplicit)
  void initSolver();

  /// initialize the stuff for the implicit solver
  void initImplicit(bool dynamic = false);

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

private:
  /// finish the computation of residual to solve in increment
  void updateResidualInternal();

  /// compute A and solve @f[ A\delta u = f_ext - f_int @f]
  template<NewmarkBeta::IntegrationSchemeCorrectorType type>
  void solveDynamic(Vector<Real> & increment);


  /* ------------------------------------------------------------------------ */
  /* Boundaries (solid_mechanics_model_boundary.cc)                           */
  /* ------------------------------------------------------------------------ */
public:
  class SurfaceLoadFunctor {
  public:
    virtual void operator()(__attribute__ ((unused)) const types::Vector<Real> & position,
			    __attribute__ ((unused)) types::Vector<Real> & force,
			    __attribute__ ((unused)) const types::Vector<Real> & normal,
			    __attribute__ ((unused)) Surface surface_id) {
      AKANTU_DEBUG_TO_IMPLEMENT();
    }

    virtual void operator()(__attribute__ ((unused)) const types::Vector<Real> & position,
			    __attribute__ ((unused)) types::Matrix & stress,
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
  /// read the material files to instantiate all the materials
  void readMaterials(const std::string & filename);

  /// read a custom material with a keyword and class as template
  template <typename M>
  UInt readCustomMaterial(const std::string & filename,
				const std::string & keyword);

  /// Use a UIntData in the mesh to specify the material to use per element
  void setMaterialIDsFromIntData(const std::string & data_name);

private:
  /// read properties part of a material file and create the material
  template <typename M>
  Material * readMaterialProperties(std::ifstream & infile,
				    ID mat_id,
				    UInt &current_line);

  /* ------------------------------------------------------------------------ */
  /* Mass (solid_mechanics_model_mass.cc)                                     */
  /* ------------------------------------------------------------------------ */
public:

  /// assemble the lumped mass matrix
  void assembleMassLumped();

  /// assemble the mass matrix
  void assembleMass();


private:
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

  inline virtual UInt getNbDataToPack(const Element & element,
				      SynchronizationTag tag) const;

  inline virtual UInt getNbDataToUnpack(const Element & element,
					SynchronizationTag tag) const;

  inline virtual void packData(CommunicationBuffer & buffer,
			       const Element & element,
			       SynchronizationTag tag) const;

  inline virtual void unpackData(CommunicationBuffer & buffer,
				 const Element & element,
				 SynchronizationTag tag);

  inline virtual UInt getNbDataToPack(SynchronizationTag tag) const;

  inline virtual UInt getNbDataToUnpack(SynchronizationTag tag) const;

  inline virtual void packData(CommunicationBuffer & buffer,
			       const UInt index,
			       SynchronizationTag tag) const;

  inline virtual void unpackData(CommunicationBuffer & buffer,
				 const UInt index,
				 SynchronizationTag tag);

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

  /**
   * @brief  get  the  SolidMechanicsModel::current_position vector  \warn  only
   * consistent after a call to SolidMechanicsModel::updateCurrentPosition
   */
  AKANTU_GET_MACRO(CurrentPosition, *current_position, const Vector<Real> &);

  /**
   * @brief get the  SolidMechanicsModel::increment vector \warn only consistent
   * if SolidMechanicsModel::setIncrementFlagOn has been called before
   */
  AKANTU_GET_MACRO(Increment,       *increment,        const Vector<Real> &);
  /// get the lumped SolidMechanicsModel::mass vector
  AKANTU_GET_MACRO(Mass,            *mass,             const Vector<Real> &);
  /// get the SolidMechanicsModel::velocity vector
  AKANTU_GET_MACRO(Velocity,        *velocity,               Vector<Real> &);

  /**
   * @brief get    the    SolidMechanicsModel::acceleration    vector,   updated    by
   * SolidMechanicsModel::updateAcceleration
   */
  AKANTU_GET_MACRO(Acceleration,   *acceleration,           Vector<Real> &);

  /// get the SolidMechnicsModel::incrementAcceleration vector
  AKANTU_GET_MACRO(IncrementAcceleration,   *increment_acceleration,           Vector<Real> &);

  /// get the SolidMechanicsModel::force vector (boundary forces)
  AKANTU_GET_MACRO(Force,           *force,                  Vector<Real> &);

  /// get the SolidMechanicsModel::residual vector, computed by SolidMechanicsModel::updateResidual
  AKANTU_GET_MACRO(Residual,        *residual,         const Vector<Real> &);

  /// get the SolidMechanicsModel::boundary vector
  AKANTU_GET_MACRO(Boundary,        *boundary,               Vector<bool> &);

  /// get the value of the SolidMechanicsModel::increment_flag
  AKANTU_GET_MACRO(IncrementFlag, increment_flag, bool);

  /// get the SolidMechanicsModel::element_material vector corresponding to a given akantu::ElementType
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementMaterial, element_material, UInt);

  /// get a particular material
  inline Material & getMaterial(UInt mat_index);
  inline const Material & getMaterial(UInt mat_index) const;

  /// give the number of materials
  inline UInt getNbMaterials() const { return materials.size(); };

  /// compute the stable time step
  Real getStableTimeStep();

  /// compute the potential energy
  Real getPotentialEnergy();
  /// compute the kinetic energy
  Real getKineticEnergy();

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

  inline FEM & getFEMBoundary(std::string name = "");

  /// get integrator
  AKANTU_GET_MACRO(Integrator, *integrator, const IntegrationScheme2ndOrder &);

private:
  /// compute the stable time step
  Real getStableTimeStep(const GhostType & ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

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

  /// materials of all elements
  ByElementTypeUInt element_material;

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

  /// Mesh
  Mesh & mesh;

  bool dynamic;

  bool implicit;
};

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "parser.hh"

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
