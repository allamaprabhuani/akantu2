/**
 * @file   structural_mechanics_model.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu May  5 15:31:11 2011
 *
 * @brief
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

#ifndef __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__
#define __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "model.hh"
#include "integrator_gauss.hh"
#include "shape_linked.hh"
#include "aka_types.hh"
#include "dumpable.hh"
#include "solver.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  class SparseMatrix;
}

__BEGIN_AKANTU__


struct StructuralMaterial {
  Real E;
  Real A;
  Real I;
  Real Iz;
  Real Iy;
  Real GJ;
};

class StructuralMechanicsModel : public Model, public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEEngineTemplate<IntegratorGauss, ShapeLinked, _ek_structural> MyFEEngineType;

  StructuralMechanicsModel(Mesh & mesh,
			   UInt spatial_dimension = _all_dimensions,
			   const ID & id = "structural_mechanics_model",
			   const MemoryID & memory_id = 0);

  virtual ~StructuralMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize fully the model
  void initFull(std::string material = "");

  /// initialize the internal vectors
  void initArrays();

  /// initialize the model
  void initModel();

  /// initialize the solver
  void initSolver(SolverOptions & options = _solver_no_options);

  /// compute the stresses per elements
  void computeStresses();

  /// assemble the stiffness matrix
  void assembleStiffnessMatrix();

  /// update the residual vector
  void updateResidual();

  /// solve the system
  void solve();


  bool testConvergenceIncrement(Real tolerance);
  bool testConvergenceIncrement(Real tolerance, Real & error);

  virtual void printself(std::ostream & stream, int indent = 0) const {};

  void computeRotationMatrix(const ElementType & type);

protected:
  UInt getTangentStiffnessVoigtSize(const ElementType & type);

  /// compute Rotation Matrices
  template<const ElementType type>
  void computeRotationMatrix(Array<Real> & rotations) {};

  /* ------------------------------------------------------------------------ */

private:
  template<ElementType type>
  inline UInt getTangentStiffnessVoigtSize();

  template <ElementType type>
  void assembleStiffnessMatrix();

  template<ElementType type>
  void computeStressOnQuad();

  template <ElementType type>
  void computeTangentModuli(Array<Real> & tangent_moduli);

  template <ElementType type>
  void transferBMatrixToSymVoigtBMatrix(Array<Real> & B, bool local = false);

  template <ElementType type>
  void transferNMatrixToSymVoigtNMatrix(Array<Real> & N_matrix);
  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);
  virtual void addDumpFieldVectorToDumper(const std::string & dumper_name,
					  const std::string & field_id);
  virtual void addDumpFieldTensorToDumper(const std::string & dumper_name,
					  const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the StructuralMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement_rotation, Array<Real> &);

  /// get the StructuralMechanicsModel::force vector (boundary forces)
  AKANTU_GET_MACRO(Force,        *force_momentum,        Array<Real> &);

 /// get the StructuralMechanicsModel::residual vector, computed by StructuralMechanicsModel::updateResidual
  AKANTU_GET_MACRO(Residual,     *residual,        const Array<Real> &);
  /// get the StructuralMechanicsModel::boundary vector
  AKANTU_GET_MACRO(BlockedDOFs,  *blocked_dofs,    Array<bool> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(RotationMatrix, rotation_matrix, Real);

  AKANTU_GET_MACRO(StiffnessMatrix, *stiffness_matrix, const SparseMatrix &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Set_ID, set_ID, UInt);

  void addMaterial(StructuralMaterial & material) { materials.push_back(material); }

  /* ------------------------------------------------------------------------ */
  /* Boundaries (structural_mechanics_model_boundary.cc)                      */
  /* ------------------------------------------------------------------------ */
public:
  /// Compute Linear load function set in global axis
  template <ElementType type>
  void computeForcesByGlobalTractionArray(const Array<Real> & tractions);

  /// Compute Linear load function set in local axis
  template <ElementType type>
  void computeForcesByLocalTractionArray(const Array<Real> & tractions);

  /// compute force vector from a function(x,y,momentum) that describe stresses
  template <ElementType type>
  void computeForcesFromFunction(BoundaryFunction in_function,
				 BoundaryFunctionType function_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// displacements array
  Array<Real> * displacement_rotation;

  /// forces array
  Array<Real> * force_momentum;

  /// stress arraz

  ElementTypeMapArray<Real> stress;

  /// residuals array
  Array<Real> * residual;

  /// boundaries array
  Array<bool> * blocked_dofs;

  /// position of a dof in the K matrix
  Array<Int> * equation_number;

  ElementTypeMapArray<UInt> element_material;

  // Define sets of beams
  ElementTypeMapArray<UInt> set_ID;

  /// local equation_number to global
  unordered_map<UInt, UInt>::type local_eq_num_to_global;

  /// stiffness matrix
  SparseMatrix * stiffness_matrix;

  /// jacobian matrix
  SparseMatrix * jacobian_matrix;

  /// increment of displacement
  Array<Real> * increment;

  /// solver for implicit
  Solver * solver;

  /// number of degre of freedom
  UInt nb_degree_of_freedom;

  // Rotation matrix
  ElementTypeMapArray<Real> rotation_matrix;


  /* -------------------------------------------------------------------------- */
  std::vector<StructuralMaterial> materials;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "structural_mechanics_model_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const StructuralMechanicsModel & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__ */
