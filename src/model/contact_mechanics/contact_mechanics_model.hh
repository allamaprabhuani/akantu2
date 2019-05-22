/**
 * @file   contact_mechanics_model.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Model of Contact Mechanics
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
#include "boundary_condition.hh"
#include "contact_detector.hh"
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_MECHANICS_MODEL_HH__
#define __AKANTU_CONTACT_MECHANICS_MODEL_HH__

namespace akantu {
class Resolution;
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class ContactMechanicsModel : public Model,
                              public DataAccessor<Element>,
                              public DataAccessor<UInt>,
                              public BoundaryCondition<ContactMechanicsModel> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

public:
  ContactMechanicsModel(
      Mesh & mesh, UInt spatial_dimension = _all_dimensions,
      const ID & id = "contact_mechanics_model", const MemoryID & memory_id = 0,
      std::shared_ptr<DOFManager> dof_manager = nullptr,
      const ModelType model_type = ModelType::_contact_mechanics_model);

  ContactMechanicsModel(
      Mesh & mesh, Array<Real> & positions,
      UInt spatial_dimension = _all_dimensions,
      const ID & id = "contact_mechanics_model", const MemoryID & memory_id = 0,
      std::shared_ptr<DOFManager> dof_manager = nullptr,
      const ModelType model_type = ModelType::_contact_mechanics_model);

  ~ContactMechanicsModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize completely the model
  void initFullImpl(const ModelOptions & options) override;

  /// allocate all vectors
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /// initialize all internal arrays for resolutions
  void initResolutions();

  /// initialize the modelType
  void initModel() override;

  /// call back for the solver, computes the force residual
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) override;

  /// callback for the solver, this assembles different matrices
  void assembleMatrix(const ID & matrix_id) override;

  /// callback for the solver, this assembles the stiffness matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const;

  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;

  /// callback for the solver, this is called at end of solve
  void afterSolveStep() override;

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /* ------------------------------------------------------------------------ */
  /* Contact Detection                                                        */
  /* ------------------------------------------------------------------------ */
public:
  void search();

  void search(Array<Real> &);

  template <Surface id> void computeNodalAreas();

  void assembleFieldsFromContactMap();

  /* ------------------------------------------------------------------------ */
  /* Contact Resolution                                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// register an empty contact resolution of a given type
  Resolution & registerNewResolution(const ID & res_name, const ID & res_type,
                                     const ID & opt_param);

protected:
  /// register a resolution in the dynamic database
  Resolution & registerNewResolution(const ParserSection & res_section);

  /// read the resolution files to instantiate all the resolutions
  void instantiateResolutions();

  /* ------------------------------------------------------------------------ */
  /* Solver Interface                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the contact stiffness matrix
  virtual void assembleStiffnessMatrix();

  /// assembles the contant internal forces
  virtual void assembleInternalForces();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  dumper::Field * createNodalFieldReal(const std::string & field_name,
                                       const std::string & group_name,
                                       bool padding_flag) override;

  dumper::Field * createNodalFieldBool(const std::string & field_name,
                                       const std::string & group_name,
                                       bool padding_flag) override;
  void dump() override;

  virtual void dump(UInt step);

  virtual void dump(Real time, UInt step);

  virtual void dump(const std::string & dumper_name);

  virtual void dump(const std::string & dumper_name, UInt step);

  virtual void dump(const std::string & dumper_name, Real time, UInt step);

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  UInt getNbData(const Array<Element> & elements,
                 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override;

  UInt getNbData(const Array<UInt> & dofs,
                 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                  const SynchronizationTag & tag) override;

protected:
  /// contact detection class
  friend class ContactDetector;

  /// contact resolution class
  friend class Resolution;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, Model::spatial_dimension, UInt);

  /// get the ContactMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement, Array<Real> &);

  /// get  the ContactMechanicsModel::increment  vector \warn  only  consistent
  /// if ContactMechanicsModel::setIncrementFlagOn has been called before
  AKANTU_GET_MACRO(Increment, *displacement_increment, Array<Real> &);

  /// get the ContactMechanics::internal_force vector (internal forces)
  AKANTU_GET_MACRO(InternalForce, *internal_force, Array<Real> &);

  /// get the ContactMechanicsModel::external_force vector (external forces)
  AKANTU_GET_MACRO(ExternalForce, *external_force, Array<Real> &);

  /// get the ContactMechanicsModel::force vector (external forces)
  Array<Real> & getForce() {
    AKANTU_DEBUG_WARNING("getForce was maintained for backward compatibility, "
                         "use getExternalForce instead");
    return *external_force;
  }

  /// get the ContactMechanics::blocked_dofs vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<Real> &);

  /// get the ContactMechanics::gaps (contact gaps)
  AKANTU_GET_MACRO(Gaps, *gaps, Array<Real> &);

  /// get the ContactMechanics::normals (normals on slave nodes)
  AKANTU_GET_MACRO(Normals, *normals, Array<Real> &);

  /// get the ContactMechanics::areas (nodal areas)
  AKANTU_GET_MACRO(NodalArea, *nodal_area, Array<Real> &);

  /// get the ContactMechanics::internal_force vector (internal forces)
  AKANTU_GET_MACRO(CurrentPositions, current_positions, Array<Real> &);

  /// get the contat map
  inline std::map<UInt, ContactElement> & getContactMap() {
    return contact_map;
  }

  ///
  inline void setPositions(Array<Real> positions) {
    detector->setPositions(positions);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// tells if the resolutions are instantiated
  bool are_resolutions_instantiated;

  /// displacements array
  Array<Real> * displacement;

  /// increment of displacement
  Array<Real> * displacement_increment{nullptr};

  /// contact forces array
  Array<Real> * internal_force{nullptr};

  /// external forces array
  Array<Real> * external_force{nullptr};

  /// boundary vector
  Array<Real> * blocked_dofs{nullptr};

  /// array to store gap between slave and master
  Array<Real> * gaps{nullptr};

  /// array to store normals from master to slave
  Array<Real> * normals{nullptr};

  /// array to store tangents on the master element
  Array<Real> * tangents{nullptr};

  /// array to store nodal areas
  Array<Real> * nodal_area{nullptr};

  /// array of current position used during update residual
  Array<Real> & current_positions;

  /// contact detection
  std::unique_ptr<ContactDetector> detector;

  /// list of contact resolutions
  std::vector<std::unique_ptr<Resolution>> resolutions;

  /// mapping between resolution name and resolution internal id
  std::map<std::string, UInt> resolutions_names_to_id;

  /// mapping between slave node its respective contact element
  std::map<UInt, ContactElement> contact_map;
};

} // namespace akantu

/* ------------------------------------------------------------------------ */
/* inline functions                                                         */
/* ------------------------------------------------------------------------ */
#include "parser.hh"
#include "resolution.hh"
/* ------------------------------------------------------------------------ */

#endif /* __AKANTU_CONTACT_MECHANICS_MODEL_HH__ */
