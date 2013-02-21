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

class StructuralMechanicsModel : public Model, public Dumpable<DumperParaview> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEMTemplate<IntegratorGauss, ShapeLinked, _ek_structural> MyFEMType;

  StructuralMechanicsModel(Mesh & mesh,
			   UInt spatial_dimension = 0,
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
  void initVectors();

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
  void computeRotationMatrix(Vector<Real> & rotations) {};

  /* ------------------------------------------------------------------------ */

private:
  template<ElementType type>
  inline UInt getTangentStiffnessVoigtSize();

  template <ElementType type>
  void assembleStiffnessMatrix();

  template<ElementType type>
  void computeStressOnQuad();

  template <ElementType type>
  void computeTangentModuli(Vector<Real> & tangent_moduli);

  template <ElementType type>
  void transferBMatrixToSymVoigtBMatrix(Vector<Real> & B, bool local = false);

  template <ElementType type>
  void transferNMatrixToSymVoigtNMatrix(Vector<Real> & N_matrix);
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

  /// get the StructuralMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement_rotation, Vector<Real> &);

  /// get the StructuralMechanicsModel::force vector (boundary forces)
  AKANTU_GET_MACRO(Force,        *force_momentum,        Vector<Real> &);

 /// get the StructuralMechanicsModel::residual vector, computed by StructuralMechanicsModel::updateResidual
  AKANTU_GET_MACRO(Residual,     *residual,        const Vector<Real> &);
  /// get the StructuralMechanicsModel::boundary vector
  AKANTU_GET_MACRO(Boundary,     *boundary,              Vector<bool> &);

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
  void computeForcesByGlobalTractionVector(const Vector<Real> & tractions);

  /// Compute Linear load function set in local axis
  template <ElementType type>
  void computeForcesByLocalTractionVector(const Vector<Real> & tractions);

  /// compute force vector from a function(x,y,momentum) that describe stresses
  template <ElementType type>
  void computeForcesFromFunction(BoundaryFunction in_function,
				 BoundaryFunctionType function_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// displacements array
  Vector<Real> * displacement_rotation;

  /// forces array
  Vector<Real> * force_momentum;

  /// stress arraz

  ByElementTypeReal stress;

  /// residuals array
  Vector<Real> * residual;

  /// boundaries array
  Vector<bool> * boundary;

  /// position of a dof in the K matrix
  Vector<Int> * equation_number;

  ByElementTypeUInt element_material;

  // Define sets of beams
  ByElementTypeUInt set_ID;

  /// local equation_number to global
  unordered_map<UInt, UInt>::type local_eq_num_to_global;

  /// stiffness matrix
  SparseMatrix * stiffness_matrix;

  /// jacobian matrix
  SparseMatrix * jacobian_matrix;

  /// increment of displacement
  Vector<Real> * increment;

  /// solver for implicit
  Solver * solver;

  /// the spatial dimension
  UInt spatial_dimension;

  /// number of degre of freedom
  UInt nb_degree_of_freedom;

  // Rotation matrix 
  ByElementTypeReal rotation_matrix;


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
