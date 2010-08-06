/**
 * @file   solid_mechanics_model.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 11:51:06 2010
 *
 * @brief  Model of Solid Mechanics
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "model.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  class Material;
  class IntegrationScheme2ndOrder;
};

__BEGIN_AKANTU__

class SolidMechanicsModel : public Model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SolidMechanicsModel(UInt spatial_dimension,
		      const ModelID & id = "solid_mechanics_model",
		      const MemoryID & memory_id = 0);

  SolidMechanicsModel(Mesh & mesh,
		      UInt spatial_dimension = 0,
		      const ModelID & id = "solid_mechanics_model",
		      const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// allocate all vectors
  void initVectors();

  /// read the material files to instantiate all the materials
  void readMaterials(const std::string & filename);

  /// initialize all internal arrays for materials
  void initMaterials();

  /// initialize the model
  void initModel();

  /// assemble the lumped mass matrix
  void assembleMass();

  /// assemble the residual for the explicit scheme
  void updateResidual();

  /// compute the acceleration from the residual
  void updateAcceleration();

  /// explicit integration predictor
  void explicitPred();

  /// explicit integration corrector
  void explicitCorr();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  AKANTU_SET_MACRO(TimeStep, time_step, Real);

  AKANTU_GET_MACRO(F_M2A, f_m2a, Real);
  AKANTU_SET_MACRO(F_M2A, f_m2a, Real);

  AKANTU_GET_MACRO(Displacement, *displacement, Vector<Real> &);
  AKANTU_GET_MACRO(Mass, *mass, Vector<Real> &);
  AKANTU_GET_MACRO(Velocity, *velocity, Vector<Real> &);
  AKANTU_GET_MACRO(Acceleration, *acceleration, Vector<Real> &);
  AKANTU_GET_MACRO(Force, *force, Vector<Real> &);
  AKANTU_GET_MACRO(Residual, *residual, Vector<Real> &);
  AKANTU_GET_MACRO(Boundary, *boundary, Vector<bool> &);

  inline Vector<Real> & getStress(ElementType type);
  inline Vector<Real> & getStrain(ElementType type);
  inline Vector<UInt> & getElementMaterial(ElementType type);

  inline Material & getMaterial(UInt mat_index);

  /// compute the stable time step
  Real getStableTimeStep();

  void setPotentialEnergyFlagOn();
  void setPotentialEnergyFlagOff();

  Real getPotentialEnergy();
  Real getKineticEnergy();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// time step
  Real time_step;

  /// conversion coefficient form force/mass to acceleration
  Real f_m2a;

  /// displacements array
  Vector<Real> * displacement;

  /// lumped mass array
  Vector<Real> * mass;

  /// velocities array
  Vector<Real> * velocity;

  /// accelerations array
  Vector<Real> * acceleration;

  /// forces array
  Vector<Real> * force;

  /// residuals array
  Vector<Real> * residual;

  /// boundaries array
  Vector<bool> * boundary;

  /// stresses arrays ordered by element types
  ByConnectivityTypeReal stress;

  /// strains arrays ordered by element types
  ByConnectivityTypeReal strain;

  /// materials of all element
  ByConnectivityTypeUInt element_material;

  /// list of used materials
  std::vector<Material *> materials;

  /// integration scheme of second order used
  IntegrationScheme2ndOrder * integrator;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "solid_mechanics_model_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SolidMechanicsModel & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_HH__ */
