

#include "aka_common.hh"
#include "constitutive_law.hh"

#ifndef AKANTU_LOCAL_LAW_HH_
#define AKANTU_LOCAL_LAW_HH_

namespace akantu{

class LocalLaw :  public ConstitutiveLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  LocalLaw(PoissonModel & model, const ID & id = "");

  virtual ~LocalLaw(){};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initConstitutiveLaw() override;

  /// constitutive law for all element of a type
  void computeFlux(ElementType el_type,
		   GhostType ghost_type = _not_ghost) override;

    /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the celerity of the fastest wave in the material
  inline Real getCelerity() const override {
    return 4.* this->epsilon;
  };

  /// compute the effective capacity
  inline Real getEffectiveCapacity() const override {
    return 1.;
  };

  
private:
  /// the permitivity
  Real epsilon;
  
  /// defines if the stiffness was computed
  bool was_stiffness_assembled;
};

} // namespace akantu

#endif /* AKANTU_LOCAL_LAW_HH_ */
