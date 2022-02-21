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

#ifndef AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH_
#define AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH_

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

protected:
  /// maximum equivalent stress
  Real norm_max_equivalent_stress;
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

  ~MaterialDamageIterative() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

  /// update internal field damage
  UInt updateDamage() override;

  UInt updateDamage(UInt quad_index, Real eq_stress, ElementType el_type,
                    GhostType ghost_type);

  /// update energies after damage has been updated
  void updateEnergiesAfterDamage(ElementType el_type) override;

  /// compute the equivalent stress on each Gauss point (i.e. the max principal
  /// stress) and normalize it by the tensile strength
  virtual void
  computeNormalizedEquivalentStress(ElementType el_type,
                                    GhostType ghost_type = _not_ghost);

  /// find max and average normalized equivalent stress
  virtual void
  findMaxNormalizedEquivalentStress(ElementType el_type,
                                    GhostType ghost_type = _not_ghost);

  inline void rotateTensor(Matrix<Real> & T,
                           const Matrix<Real> & rotation_matrix);

protected:
  inline UInt updateDamageOnQuad(UInt quad_index, Real /*eq_stress*/,
                                 ElementType el_type, GhostType ghost_type);

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

  /// internal field to store equivalent stress on each Gauss point
  InternalField<Real> equivalent_stress;

  /// the reduction
  InternalField<UInt> reduction_step;

  /// 1st vector normal to crack, 2nd (& 3rd) - vector(s) in crack plane
  InternalField<Real> crack_normals;

  // array to store previously converged damage
  InternalField<Real> damage_stored;

  // array to store previously converged number of reduction steps
  InternalField<UInt> reduction_step_stored;

  /// additional volume of damaged elements
  InternalField<Real> extra_volume;

  /// damage increment
  Real prescribed_dam;

  /// define damage threshold at which damage will be set to 1
  Real dam_threshold;

  /// deviation from max stress at which Gauss point will still get damaged
  Real dam_tolerance;

  /// maximum damage value
  Real max_damage;

  /// the number of total reductions steps until complete failure
  UInt max_reductions;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_damage_iterative_inline_impl.hh"

#endif /* AKANTU_MATERIAL_DAMAGE_ITERATIVE_HH_ */
