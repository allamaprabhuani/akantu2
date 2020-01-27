/**
 * @file   material_FE2.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief Material for multi-scale simulations. It stores an
 * underlying RVE on each integration point of the material.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "material_thermal.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_FE_2_HH__
#define __AKANTU_MATERIAL_FE_2_HH__

namespace akantu {
class SolidMechanicsModelRVE;
}

namespace akantu {

/* -------------------------------------------------------------------------- */
/// /!\ This material works ONLY for meshes with a single element type!!!!!
/* -------------------------------------------------------------------------- */

/**
 * MaterialFE2
 *
 * parameters in the material files :
 *   - mesh_file
 */
template <UInt DIM> class MaterialFE2 : public MaterialThermal<DIM> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef MaterialThermal<DIM> Parent;

public:
  MaterialFE2(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialFE2();

  typedef VoigtHelper<DIM> voigt_h;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initMaterial();

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

  /// gel strain inrease linear with time, exponential with temperature
  void computeNewGelStrain(Matrix<Real> & gelstrain, const Real & delta_time,
                           const Real & temp);

  /// assymptotic gel strain - time curve
  void computeNewGelStrainTimeDependent(Matrix<Real> & gelstrain,
                                        const Real & delta_time, const Real & T,
                                        Real & non_reacted_gel);

  /// advance alkali-silica reaction by the user-provided gel strain
  void advanceASR(const Matrix<Real> & prestrain);

  /// advance alkali-silica reaction based on delta time and temperature-
  /// dependent reaction rate
  void advanceASR(const Real & delta_time);

  /// compute amount of gel strain averaged across all RVEs
  Real computeAverageGelStrain();

  /// set default dumper directory to all rves
  void setDirectoryToRveDumper(const std::string & directory);

  /// compute number of RVEs
  UInt getNbRVEs();

  /// dump all the rves
  void dump();

  /// increase gel strain according to time step
  void increaseGelStrain(Real & dt);

  /// update damage ratio after converged step
  virtual void afterSolveStep();

private:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Underlying RVE at each integration point
  std::vector<std::unique_ptr<SolidMechanicsModelRVE>> RVEs;

  /// Meshes for all RVEs
  std::vector<std::unique_ptr<Mesh>> meshes;

  /// the element type of the associated mesh (this material handles only one
  /// type!!)
  ElementType el_type;

  /// the name of RVE mesh file
  ID mesh_file;

  /// Elastic stiffness tensor at each Gauss point (in voigt notation)
  InternalField<Real> C;

  /// number of gel pockets in each underlying RVE
  UInt nb_gel_pockets;

  /// pre-exponential factor of Arrhenius law
  Real k;

  /// activation energy of ASR in Arrhenius law
  Real activ_energy;

  /// universal gas constant;
  Real R;

  /// saturation constant for time dependent gel strain increase
  Real sat_const;

  /// current gelstrain due to ASR at each Gauss point
  InternalField<Real> gelstrain;

  /// percent of yet non-reacted gel (for time-dependent asr simulation)
  InternalField<Real> non_reacted_gel;

  /// ratio between area of damaged elements weighted by damage value
  /// and the total area of RVE
  InternalField<Real> damage_ratio;

  InternalField<Real> damage_ratio_paste;
  InternalField<Real> damage_ratio_agg;

  /// Macro eigen stress incrementally summed from homogenized meso-scale
  InternalField<Real> eigen_stress;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_FE2_inline_impl.cc"

} // namespace akantu

#endif /* __AKANTU_MATERIAL_FE_2_HH__ */
