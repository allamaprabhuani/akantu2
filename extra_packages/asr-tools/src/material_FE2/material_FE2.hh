
/**
 * @file   material_FE2.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief Material for multi-scale simulations. It stores an
 * underlying RVE on each integration point of the material.
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "material_thermal.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_FE_2_HH_
#define AKANTU_MATERIAL_FE_2_HH_

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
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

  // /// compute ASR strain according to the sigmoidal rule (Larive,1998)
  // void computeASRStrainLarive(const Real & delta_time_day, const Real & T,
  //                             Matrix<Real> & gelstrain);
  // /// gel strain inrease linear with time, exponential with temperature
  // void computeNewGelStrain(Matrix<Real> & gelstrain,
  //                          const Real & delta_time_day, const Real &
  //                          temp);

  // /// assymptotic gel strain - time curve
  // void computeNewGelStrainTimeDependent(Matrix<Real> & gelstrain,
  //                                       const Real & delta_time_day,
  //                                       const Real & T, Real &
  //                                       non_reacted_gel);

  // /// reset the gel strain value to a previous value
  // void resetGelStrain(const Real & old_time_step);

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

  /// dump all the rves with the dumping number
  void dump(UInt dump_nb);

  /// increase gel strain according to time step
  void increaseGelStrain(Real & dt);

  /// set time step to all RVEs
  void setTimeStep(Real dt);

  /// update damage ratio after converged step
  virtual void afterSolveStep();

  /// save state of all rves
  void saveRVEsState(std::string & output_dir);

  /// load state of all rves
  void loadRVEsState(std::string & output_dir);

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

  /// parameters of sigmoidal expansion law (Larive, 1998)
  Real eps_inf, U_C, U_L, T_ref, time_lat_ref, time_ch_ref;

  /// current gelstrain due to ASR at each Gauss point
  InternalField<Real> gelstrain;

  /// ratio between area of damaged elements weighted by damage value
  /// and the materials' areas
  InternalField<Real> crack_volume_ratio;
  InternalField<Real> crack_volume_ratio_paste;
  InternalField<Real> crack_volume_ratio_agg;

  /// flag to reset damage to previously converged values on each iteration
  bool reset_damage{false};
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_MATERIAL_FE_2_HH_ */
