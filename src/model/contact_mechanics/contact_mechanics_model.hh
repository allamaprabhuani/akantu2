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
#include "data_accessor.hh"
#include "model.hh"
#include "fe_engine.hh"
#include "contact_detector.hh"
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
class ContactMechanicsModel :
    public Model,
    public DataAccessor<Element>,
    public DataAccessor<UInt> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
 
  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;
  
public:
  ContactMechanicsModel(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
			const ID & id = "contact_mechanics_model",
			const MemoryID & memory_id = 0,
			const ModelType model_type = ModelType::_contact_mechanics_model);

  ContactMechanicsModel(Mesh & mesh, Array<Real> & positions, UInt spatial_dimension = _all_dimensions,
			const ID & id = "contact_mechanics_model",
			const MemoryID & memory_id = 0,
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
  
  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /* ------------------------------------------------------------------------ */
  /* Contact Detection                                                        */
  /* ------------------------------------------------------------------------ */
public:
  void search();
  
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
protected:
  
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
  
  virtual void dump(const std::string & dumper_name);

  virtual void dump(const std::string & dumper_name, UInt step);

  virtual void dump(const std::string & dumper_name, Real time, UInt step);

  void dump() override;

  virtual void dump(UInt step);

  virtual void dump(Real time, UInt step);

  
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
  /// get the ContactMechanics::contact_force vector (internal forces)
  AKANTU_GET_MACRO(InternalForce, *contact_force, Array<Real> &);

  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<Real> &);
 
  AKANTU_GET_MACRO(Gaps, *gaps, Array<Real> &);

  /// get the contat map
  inline std::map<UInt, ContactElement> & getContactMap() { return contact_map; }
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// tells if the resolutions are instantiated
  bool are_resolutions_instantiated;

  /// contact forces array
  Array<Real> * contact_force{nullptr};

  /// boundary vector
  Array<Real> * blocked_dofs{nullptr};

  /// gaps vector
  Array<Real> * gaps{nullptr};
 
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
#include "resolution.hh"
#include "parser.hh"
/* ------------------------------------------------------------------------ */



#endif  /* __AKANTU_CONTACT_MECHANICS_MODEL_HH__ */
