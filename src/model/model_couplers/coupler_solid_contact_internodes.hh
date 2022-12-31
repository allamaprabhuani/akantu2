/* -------------------------------------------------------------------------- */
#include "coupler_solid_contact.hh"
#include "contact_mechanics_internodes_model.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_COUPLER_SOLID_CONTACT_INTERNODES_HH__
#define __AKANTU_COUPLER_SOLID_CONTACT_INTERNODES_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
class CouplerSolidContactInternodesTemplate :
    public AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsInternodesModel> {
public:
  CouplerSolidContactInternodesTemplate(
      Mesh & mesh, UInt dim = _all_dimensions,
      const ID & id = "coupler_solid_contact_internodes",
      std::shared_ptr<DOFManager> dof_manager = nullptr);

  ~CouplerSolidContactInternodesTemplate() override = default;

  /// custom solve step to handle the INTERNODES iterations
  void solveStep(const ID & solver_id = "") override;

protected:
  /// callback for the solver to assemble the rhs
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) const override;

  /// callback for the solver, this assembles different matrices
  void assembleMatrix(const ID & matrix_id) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;
};

/* -------------------------------------------------------------------------- */
using CouplerSolidContactInternodes =
    CouplerSolidContactInternodesTemplate<SolidMechanicsModel>;

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "coupler_solid_contact_internodes_tmpl.hh"

#endif // __AKANTU_COUPLER_SOLID_CONTACT_INTERNODES_HH__
