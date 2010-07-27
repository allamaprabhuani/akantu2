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
#include "common.hh"
#include "model.hh"
#include "material.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class SolidMechanicsModel : public Model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SolidMechanicsModel(UInt spatial_dimension,
	const ModelID & id = "solid_mechanics_model",
	const MemoryID & memory_id = 0);

  SolidMechanicsModel(UInt spatial_dimension,
	Mesh & mesh,
	const ModelID & id = "solid_mechanics_model",
	const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize all internal arrays, mesh must be set before
  void initModel();

  /// assemble the lumped mass matrix
  void assembleMass();

  /// read the material files to instantiate all the materials
  void readMaterials(const std::string & filename);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// displacements arrays ordered by element types
  ByConnectivityTypeReal displacement;

  /// lumped mass arrays ordered by element types
  ByConnectivityTypeReal mass;

  /// velocities arrays ordered by element types
  ByConnectivityTypeReal velocity;

  /// accelerations arrays ordered by element types
  ByConnectivityTypeReal acceleration;

  /// forces arrays ordered by element types
  ByConnectivityTypeReal force;

  /// residuals arrays ordered by element types
  ByConnectivityTypeReal residual;

  /// stresses arrays ordered by element types
  ByConnectivityTypeReal stress;

  /// strains arrays ordered by element types
  ByConnectivityTypeReal strain;

  /// boundaries arrays ordered by element types
  ByConnectivityTypeInt boundary;

  /// materials of all element
  ByConnectivityTypeInt element_material;

  /// list of used materials
  std::vector<MaterialBase> materials;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

// #include "solid_mechanics_model_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SolidMechanicsModel & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_HH__ */
