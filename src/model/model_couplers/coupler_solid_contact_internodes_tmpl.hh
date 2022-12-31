/* -------------------------------------------------------------------------- */
#include "coupler_solid_contact_internodes.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactInternodesTemplate<SolidMechanicsModelType>
    ::solveStep(const ID & solver_id) {
  const UInt MAX_STEPS = 100; // security value for now, increase later if needed

  UInt step = 0;
  while (true) {
    // Always search before every step
    this->contact->search();

    // Perform regular solve
    Model::solveStep(solver_id);

    if (++step >= MAX_STEPS) {
      AKANTU_EXCEPTION("Reached maximum number of internodes steps: " << step);
    }

    // Perform any required update, and solve again if needed
    auto & current_positions = this->contact->getContactDetector().getPositions();
    current_positions.copy(this->solid->getCurrentPosition());
    if (!this->contact->updateAfterStep()) {
      break;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
MatrixType CouplerSolidContactInternodesTemplate<SolidMechanicsModelType>
    ::getMatrixType(const ID & matrix_id) const {
  if (matrix_id == "K") {
    return _unsymmetric;
  }
  // TODO: this wasn't in the original internodes model, but it should be symm?
//  if (matrix_id == "M") {
//    return _symmetric;
//  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactInternodesTemplate<SolidMechanicsModelType>
    ::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    // We need the solid's mass matrix to compute the stiffness matrix.
    this->solid->assembleMass();
    this->contact->assembleInternodesMatrix();
  } else if (matrix_id == "M") {
    this->solid->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
ModelSolverOptions
CouplerSolidContactInternodesTemplate<SolidMechanicsModelType>
    ::getDefaultSolverOptions(const TimeStepSolverType & type) const {
  return this->contact->getDefaultSolverOptions(type);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactInternodesTemplate<SolidMechanicsModelType>
    ::assembleResidual() {
//  switch (this->method) {
  // TODO: this is needed because contact search hasn't been performed before,
  // but why does the ContactMechanicsModel coupler not do it then? (_static)
//  case _static:
//  case _explicit_lumped_mass: {
//    auto & current_positions = this->contact->getContactDetector().getPositions();
//    current_positions.copy(this->solid->getCurrentPosition());
//    this->contact->search();
//    break;
//  }
//  default:
//    break;
//  }

  this->solid->assembleResidual();
  this->contact->assembleResidual();
}

} // namespace akantu