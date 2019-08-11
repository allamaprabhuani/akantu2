/**
 * @file   material_damage_iterative.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Specialization of the class material damage to damage only one gauss
 * point at a time and propagate damage in a linear way. Max principal stress
 * criterion is used as a failure criterion.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH__
#define __AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH__

namespace akantu {

class MaterialDamageIterativeInterface {
public:
  virtual ~MaterialDamageIterativeInterface() = default;
  virtual UInt updateDamage() = 0;
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get max normalized equivalent stress
  AKANTU_GET_MACRO(NormMaxEquivalentStress, norm_max_equivalent_stress, Real);

  /// get a non-const reference to the max normalized equivalent stress
  AKANTU_GET_MACRO_NOT_CONST(NormMaxEquivalentStress,
                             norm_max_equivalent_stress, Real &);

  /// get average normalized equivalent stress
  AKANTU_GET_MACRO(NormAvEquivalentStress, norm_av_equivalent_stress, Real);

protected:
  /// maximum equivalent stress
  Real norm_max_equivalent_stress;

  /// average equivalent stress
  Real norm_av_equivalent_stress;
};

/**
 * Material damage iterative
 *
 * parameters in the material files :
 *   - Sc
 */
template <UInt spatial_dimension,
          template <UInt> class ElasticParent = MaterialElastic>
class MaterialDamageIterative
    : public MaterialDamage<spatial_dimension, ElasticParent>,
      public MaterialDamageIterativeInterface {
  using parent = MaterialDamage<spatial_dimension, ElasticParent>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDamageIterative(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialDamageIterative(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

  /// update internal field damage
  UInt updateDamage() override;

  UInt updateDamageOnQuad(UInt quad_index, const Real eq_stress,
                          const ElementType & el_type,
                          const GhostType & ghost_type);

  /// update energies after damage has been updated
  virtual void updateEnergiesAfterDamage(ElementType el_type);

  /// compute the equivalent stress on each Gauss point (i.e. the max prinicpal
  /// stress) and normalize it by the tensile strength
  virtual void
  computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                    ElementType el_type,
                                    GhostType ghost_type = _not_ghost);

  /// find max and average normalized equivalent stress
  virtual void
  findMaxAndAvNormalizedEquivalentStress(ElementType el_type,
                                         GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma, Real & dam);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */

  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// resistance to damage
  RandomInternalField<Real> Sc;

  /// the reduction
  InternalField<UInt> reduction_step;

  /// internal field to store equivalent stress on each Gauss point
  InternalField<Real> equivalent_stress;

  /// the number of total reductions steps until complete failure
  UInt max_reductions;

  /// damage increment
  Real prescribed_dam;

  /// deviation from max stress at which Gauss point will still get damaged
  Real dam_tolerance;

  /// define damage threshold at which damage will be set to 1
  Real dam_threshold;

  /// maximum damage value
  Real max_damage;

  /// recovery of stiffness when in compression
  bool contact;

  /// trace of a stress tensor to judge on compression
  InternalField<Real> stress_trace;

};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_damage_iterative_inline_impl.cc"

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH__ */
