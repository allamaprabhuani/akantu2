/* -------------------------------------------------------------------------- */
#include "phasefield.hh"
#include "volumetric_deviatoric_split.hh"
#include "no_energy_split.hh"
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_LINEAR_HH__
#define __AKANTU_PHASEFIELD_LINEAR_HH__

namespace akantu {

template <Int dim> class PhaseFieldLinear : public PhaseField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PhaseFieldLinear(PhaseFieldModel & model, const ID & id = "");

  ~PhaseFieldLinear() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the dissiapted energy
  void computeDissipatedEnergy(ElementType el_type) override;

  void
  computeDissipatedEnergyByElement(const Element & element,
                                   Vector<Real> & edis_on_quad_points) override;

  void afterSolveStep() override;

protected:
  void
  computeDissipatedEnergyByElement(ElementType type, Idx index,
                                   Vector<Real> & edis_on_quad_points) override;

  void computePhiOnQuad(const Matrix<Real> & /*strain_quad*/,
                        Real & /*phi_quad*/);

  void computePhiIsotropicOnQuad(const Matrix<Real> & /*strain_quad*/,
                                 Real & /*phi_quad*/);

  void computeDrivingForce(ElementType /*el_type*/,
                           GhostType /*ghost_type*/) override;

  inline void computeDrivingForceOnQuad(const Real & /*phi_quad*/,
                                        Real & /*driving_force_quad*/,
                                        const Real & /*g_c_quad*/);

  inline void computeDamageEnergyDensityOnQuad(const Real & /*phi_quad*/,
                                               Real & /*dam_energy_quad*/);

  inline void
  computeDissipatedEnergyOnQuad(const Real & /*dam_quad*/,
                                const Real & /*dam_prev_quad*/,
                                const Vector<Real> & /*grad_d_quad */,
                                Real & /*energy*/, Real & /*g_c_quad*/);

public:
  void updateInternalParameters() override;

  void initPhaseField() override;

private:
  // irreversibility tolerance
  Real tol_ir;

  // penalization parameter
  Real gamma;

  // dimension to consider in deviatoric split
  Int dev_dim;

  // energy split
  std::shared_ptr<EnergySplit> energy_split{nullptr};
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "phasefield_linear_inline_impl.hh"

#endif
