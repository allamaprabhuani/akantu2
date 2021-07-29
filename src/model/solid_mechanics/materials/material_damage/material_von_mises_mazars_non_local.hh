


/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_damage_non_local.hh"
#include "material_von_mises_mazars.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_VON_MISES_MAZARS_NON_LOCAL_HH_
#define AKANTU_MATERIAL_VON_MISES_MAZARS_NON_LOCAL_HH_

namespace akantu {

/**
 * Material Mazars Non local + Von Mises plasticity
 *
 * parameters in the material files :
 */
template <UInt spatial_dimension>
class MaterialVonMisesMazarsNonLocal
    : public MaterialDamageNonLocal<spatial_dimension,
                                    MaterialVonMisesMazars<spatial_dimension>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MaterialNonLocalParent =
      MaterialDamageNonLocal<spatial_dimension,
                             MaterialVonMisesMazars<spatial_dimension>>;

  MaterialVonMisesMazarsNonLocal(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  void computeNonLocalStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost) override;

  void registerNonLocalVariables() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the ehat per quadrature points to perform the averaging
  InternalField<Real> Ehat;

  InternalField<Real> non_local_variable;
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_VON_MISES_MAZARS_NON_LOCAL_HH_ */
