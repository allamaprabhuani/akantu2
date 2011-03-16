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
#include "material.hh"
#include "material_parser.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "aka_types.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  //  class Material;
  class IntegrationScheme2ndOrder;
  class Contact;
  class Solver;
  class SparseMatrix;
}

__BEGIN_AKANTU__

class SolidMechanicsModel : public Model {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEMTemplate<IntegratorGauss,ShapeLagrange> MyFEMType;

  SolidMechanicsModel(UInt spatial_dimension,
		      const ModelID & id = "solid_mechanics_model",
		      const MemoryID & memory_id = 0);

  SolidMechanicsModel(Mesh & mesh,
		      UInt spatial_dimension = 0,
		      const ModelID & id = "solid_mechanics_model",
		      const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// allocate all vectors
  void initVectors();

  /// initialize all internal arrays for materials
  void initMaterials();

  /// initialize the model
  void initModel();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Explicit                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the array needed by updateResidual (residual, current_position)
  void initializeUpdateResidualData();

  /// update the current position vector
  void updateCurrentPosition();

  /// assemble the residual for the explicit scheme
  void updateResidual(bool need_initialize = true);

  /// compute the acceleration from the residual
  void updateAcceleration();

  /// explicit integration predictor
  void explicitPred();

  /// explicit integration corrector
  void explicitCorr();

  /* ------------------------------------------------------------------------ */
  /* Implicit                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the stuff for the implicit solver
  void initImplicitSolver();

  /// assemble the stiffness matrix
  void assembleStiffnessMatrix();

  /// solve Ku = f
  //  void solve(Vector<Real> & solution);
  void solve();

  /// test the convergence (norm of increment)
  bool testConvergenceIncrement(Real tolerance);

  /// test the convergence (norm of residual)
  bool testConvergenceResidual(Real tolerance);

  /* ------------------------------------------------------------------------ */
  /* Boundaries (solid_mechanics_model_boundary.cc)                           */
  /* ------------------------------------------------------------------------ */
public:
  /// integrate a force on the boundary by providing a stress tensor
  void computeForcesByStressTensor(const Vector<Real> & stresses,
				   const ElementType & type);

  /// integrate a force on the boundary by providing a traction vector
  void computeForcesByTractionVector(const Vector<Real> & tractions,
				     const ElementType & type);

  /// compute force vector from a function(x,y,z) that describe stresses
  void computeForcesFromFunction(BoundaryFunction in_function,
				 BoundaryFunctionType function_type);

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

private:
  /// read properties part of a material file and create the material
  template <typename M>
  Material * readMaterialProperties(std::ifstream & infile,
				    MaterialID mat_id,
				    UInt &current_line);

  /* ------------------------------------------------------------------------ */
  /* Mass (solid_mechanics_model_mass.cc)                                     */
  /* ------------------------------------------------------------------------ */
public:

  /// assemble the lumped mass matrix
  void assembleMassLumped();

private:
  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumped(GhostType ghost_type);

  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumpedRowSum(GhostType ghost_type, const ElementType type);

  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumpedDiagonalScaling(GhostType ghost_type, const ElementType type);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
public:

  inline virtual UInt getNbDataToPack(const Element & element,
				      GhostSynchronizationTag tag) const;

  inline virtual UInt getNbDataToUnpack(const Element & element,
					GhostSynchronizationTag tag) const;

  inline virtual void packData(Real ** buffer,
			       const Element & element,
			       GhostSynchronizationTag tag) const;

  inline virtual void unpackData(Real ** buffer,
				 const Element & element,
				 GhostSynchronizationTag tag) const;

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
  /// get the SolidMechanicsModel::force vector (boundary forces)
  AKANTU_GET_MACRO(Force,           *force,                  Vector<Real> &);
  /// get the SolidMechanicsModel::residual vector, computed by SolidMechanicsModel::updateResidual
  AKANTU_GET_MACRO(Residual,        *residual,         const Vector<Real> &);
  /// get the SolidMechanicsModel::boundary vector
  AKANTU_GET_MACRO(Boundary,        *boundary,               Vector<bool> &);

  /// get the value of the SolidMechanicsModel::increment_flag
  AKANTU_GET_MACRO(IncrementFlag, increment_flag, bool);

  /// get the SolidMechanicsModel::element_material vector corresponding to a given akantu::ElementType
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, Vector<UInt> &);

  /// get a particular material
  inline Material & getMaterial(UInt mat_index);

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

  /// get the equation number Vector<Int>
  AKANTU_GET_MACRO(EquationNumber, *equation_number, const Vector<Int> &);

  /// get the stiffness matrix
  AKANTU_GET_MACRO(StiffnessMatrix, *stiffness_matrix, SparseMatrix &);

  inline FEM & getFEMBoundary(std::string name = "");

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

  /// forces array
  Vector<Real> * force;

  /// residuals array
  Vector<Real> * residual;

  /// boundaries array
  Vector<bool> * boundary;

  /// array of current position used during update residual
  Vector<Real> * current_position;

  /// position of a dof in the K matrix
  Vector<Int> * equation_number;

  /// local equation_number to global
  unordered_map<UInt, UInt>::type local_eq_num_to_global;

  /// stiffness matrix
  SparseMatrix * stiffness_matrix;

  /// materials of all elements
  ByElementTypeUInt element_material;

  /// materials of all ghost elements
  ByElementTypeUInt ghost_element_material;

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

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "solid_mechanics_model_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SolidMechanicsModel & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_HH__ */
