

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_linear_isotropic_hardening.hh"
#include "material_mazars.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VONMISES_MAZARS_HH__
#define __AKANTU_MATERIAL_VONMISES_MAZARS_HH__

namespace akantu {

template <UInt spatial_dimension>
class MaterialVonMisesMazars
  : public MaterialLinearIsotropicHardening<spatial_dimension>,
    public MaterialMazars<spatial_dimension> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialVonMisesMazars(SolidMechanicsModel & model,
			 const ID & id = "");
  MaterialVonMisesMazars(SolidMechanicsModel & model, UInt dim,
			 const Mesh & mesh, FEEngine & fe_engine,
			 const ID & id = "");
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override; 

};

}

#endif /* __AKANTU_MATERIAL_VONMISES_MAZARS_HH__ */

