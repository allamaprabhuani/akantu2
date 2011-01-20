/**
 * @file   material.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jul 23 09:06:29 2010
 *
 * @brief  Mother class for all materials
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_HH__
#define __AKANTU_MATERIAL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "fem.hh"
#include "mesh.hh"
//#include "solid_mechanics_model.hh"
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  class SolidMechanicsModel;
}

__BEGIN_AKANTU__

class Material : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Material(SolidMechanicsModel & model, const MaterialID & id = "");
  virtual ~Material();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// read properties
  virtual void setParam(const std::string & key, const std::string & value,
			const MaterialID & id);

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// compute the residual for this material
  void updateResidual(Vector<Real> & current_position,
		      GhostType ghost_type = _not_ghost);

  /// constitutive law
  virtual void computeStress(ElementType el_type,
			     GhostType ghost_type = _not_ghost) = 0;

  /// compute the potential energy by element
  void computePotentialEnergyByElement();

  /// compute the potential energy
  virtual void computePotentialEnergy(ElementType el_type,
				      GhostType ghost_type = _not_ghost) = 0;

  /// compute the stable time step for an element of size h
  virtual Real getStableTimeStep(Real h) = 0;

  /// add an element to the local mesh filter
  inline void addElement(ElementType type, UInt element);

  /// add an element to the local mesh filter for ghost element
  inline void addGhostElement(ElementType type, UInt element);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Function for all materials                                               */
  /* ------------------------------------------------------------------------ */
protected:
  /// allocate an internal vector
  void initInternalVector(ByElementTypeReal & vect,
			  UInt nb_component,
			  const std::string & id,
			  GhostType ghost_type = _not_ghost);

  /// resize an internal vector
  void resizeInternalVector(ByElementTypeReal & vect,
			    GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
public:

  inline virtual UInt getNbDataToPack(const Element & element,
				      GhostSynchronizationTag tag);

  inline virtual UInt getNbDataToUnpack(const Element & element,
					GhostSynchronizationTag tag);

  inline virtual void packData(Real ** buffer,
			       const Element & element,
			       GhostSynchronizationTag tag);

  inline virtual void unpackData(Real ** buffer,
				 const Element & element,
				 GhostSynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(ID, id, const MaterialID &);
  AKANTU_GET_MACRO(Rho, rho, Real);

  Real getPotentialEnergy();

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Strain, strain, const Vector<Real> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Stress, stress, const Vector<Real> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(PotentialEnergy, potential_energy, const Vector<Real> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the material
  MaterialID id;

  /// spatial dimension
  UInt spatial_dimension;

  /// material name
  std::string name;

  /// The model to witch the material belong
  SolidMechanicsModel * model;

  /// density : rho
  Real rho;

  /// stresses arrays ordered by element types
  ByElementTypeReal stress;

  /// strains arrays ordered by element types
  ByElementTypeReal strain;

  /// list of element handled by the material
  ByElementTypeUInt element_filter;

  /// stresses arrays ordered by element types
  ByElementTypeReal ghost_stress;

  /// strains arrays ordered by element types
  ByElementTypeReal ghost_strain;

  /// list of element handled by the material
  ByElementTypeUInt ghost_element_filter;

  /// is the vector for potential energy initialized
  bool potential_energy_vector;

  /// potential energy by element
  ByElementTypeReal potential_energy;

  /// potential energy by element
  ByElementTypeReal ghost_potential_energy;

  /// boolean to know if the material has been initialized
  bool is_init;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Material & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#include "materials/material_elastic.hh"

/* -------------------------------------------------------------------------- */
/* Auto loop                                                                  */
/* -------------------------------------------------------------------------- */

#define MATERIAL_QUADRATURE_POINT_LOOP_BEGIN				\
  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(el_type);	\
  UInt size_strain          = spatial_dimension * spatial_dimension;	\
  									\
  UInt nb_element;							\
  Real * strain_val;							\
  Real * stress_val;							\
  									\
  if(ghost_type == _not_ghost) {					\
    nb_element   = element_filter[el_type]->getSize();			\
    stress[el_type]->resize(nb_element*nb_quadrature_points);		\
    strain_val = strain[el_type]->values;				\
    stress_val = stress[el_type]->values;				\
  } else {								\
    nb_element = ghost_element_filter[el_type]->getSize();		\
    ghost_stress[el_type]->resize(nb_element*nb_quadrature_points);	\
    strain_val = ghost_strain[el_type]->values;				\
    stress_val = ghost_stress[el_type]->values;				\
  }									\
  									\
  if (nb_element == 0) return;						\
  									\
  for (UInt el = 0; el < nb_element; ++el) {				\
    for (UInt q = 0; q < nb_quadrature_points; ++q) {			\


#define MATERIAL_QUADRATURE_POINT_LOOP_END				\
      strain_val += size_strain;					\
      stress_val += size_strain;					\
    }									\
  }                                                                     \


#endif /* __AKANTU_MATERIAL_HH__ */
