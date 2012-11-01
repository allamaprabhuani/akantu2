/**
 * @file   structural_mechanics_model.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu May  5 15:31:11 2011
 *
 * @brief StructuralMechanicsModel description
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

/* -------------------------------------------------------------------------- */
namespace akantu {
  class Solver;
  class SparseMatrix;
}

__BEGIN_AKANTU__


struct StructuralMaterial {
  Real E;
  Real A;
  Real I;
};

class StructuralMechanicsModel : public Model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEMTemplate<IntegratorGauss,ShapeLinked> MyFEMType;

  StructuralMechanicsModel(UInt spatial_dimension,
		      const ID & id = "structural_mechanics_model",
		      const MemoryID & memory_id = 0);

  StructuralMechanicsModel(Mesh & mesh,
			   UInt spatial_dimension = 0,
			   const ID & id = "structural_mechanics_model",
			   const MemoryID & memory_id = 0);

  virtual ~StructuralMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initVectors();
  void initModel();
  void initImplicitSolver();
  void computeStressOnQuad();
  void assembleStiffnessMatrix();
  void updateResidual();
  void solve();
  bool testConvergenceIncrement(Real tolerance);
  bool testConvergenceIncrement(Real tolerance, Real & error);
  virtual void printself(__attribute__ ((unused)) std::ostream & stream,
			 __attribute__ ((unused)) int indent = 0) const {};
  UInt getTangentStiffnessVoigtSize(const ElementType & type);

  /* ------------------------------------------------------------------------ */

private:
  template<ElementType type>
  inline UInt getTangentStiffnessVoigtSize();

  template <ElementType type>
  void assembleStiffnessMatrix();

  template<ElementType type>
  void computeStressOnQuad();

  template <ElementType type>
  void computeTangentStiffness(Vector<Real> & tangent_stiffness_matrix);

  template <ElementType type>
  void transferBMatrixToSymVoigtBMatrix(Vector<Real> & B);

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

  AKANTU_GET_MACRO(StiffnessMatrix, *stiffness_matrix, const SparseMatrix &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt)

  void addMaterial(StructuralMaterial & material) { materials.push_back(material); }

  /* ------------------------------------------------------------------------ */
  /* Boundaries (structural_mechanics_model_boundary.cc)                      */
  /* ------------------------------------------------------------------------ */
public:
  /// !!NOT IMPLEMENTED YET!!
  void computeForcesByStressTensor(const Vector<Real> & stresses,
				   const ElementType & type);

  /// integrate a force on the boundary by providing a traction vector
  void computeForcesByTractionVector(const Vector<Real> & tractions,
				     const ElementType & type);

  /// compute force vector from a function(x,y,momentum) that describe stresses
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

  /// local equation_number to global
  unordered_map<UInt, UInt>::type local_eq_num_to_global;

  /// stiffness matrix
  SparseMatrix * stiffness_matrix;

  /// increment of displacement
  Vector<Real> * increment;

  /// solver for implicit
  Solver * solver;

  /// the spatial dimension
  UInt spatial_dimension;

  /// number of degre of freedom
  UInt nb_degree_of_freedom;

  /* -------------------------------------------------------------------------- */
  std::vector<StructuralMaterial> materials;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "structural_mechanics_model_inline_impl.cc"
#endif

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const StructuralMechanicsModel & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__ */
