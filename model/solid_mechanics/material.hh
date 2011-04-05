/**
 * @file   material.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jul 23 09:06:29 2010
 *
 * @brief  Mother class for all materials
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
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


  /// compute the stiffness matrix
  void assembleStiffnessMatrix(Vector<Real> & current_position,
			       GhostType ghost_type);

  /// compute the stable time step for an element of size h
  virtual Real getStableTimeStep(Real h) = 0;

  /// add an element to the local mesh filter
  inline void addElement(ElementType type, UInt element);

  /// add an element to the local mesh filter for ghost element
  inline void addGhostElement(ElementType type, UInt element);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

protected:

  /// constitutive law
  virtual void computeStress(ElementType el_type,
			     GhostType ghost_type = _not_ghost) = 0;

  /// compute the tangent stiffness matrix
  virtual void computeTangentStiffness(const ElementType & el_type,
				       Vector<Real> & tangent_matrix,
				       GhostType ghost_type = _not_ghost) = 0;

  /// compute the potential energy
  virtual void computePotentialEnergy(ElementType el_type,
				      GhostType ghost_type = _not_ghost) = 0;

private:
  template<UInt dim>
  void assembleStiffnessMatrix(Vector<Real> & current_position,
			       const ElementType & type,
			       GhostType ghost_type);

  /// transfer the B matrix to a Voigt notation B matrix
  template<UInt dim>
  inline void transferBMatrixToSymVoigtBMatrix(Real * B, Real * Bvoigt, UInt nb_nodes_per_element) const;

  inline UInt getTangentStiffnessVoigtSize(UInt spatial_dimension) const;


  /// compute the potential energy by element
  void computePotentialEnergyByElement();

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

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementFilter, element_filter, const Vector<UInt> &);
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
#include "materials/material_damage.hh"
//#include "materials/material_mazars.hh"

/* -------------------------------------------------------------------------- */
/* Auto loop                                                                  */
/* -------------------------------------------------------------------------- */

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN			\
  UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(el_type); \
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


#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END			\
      strain_val += size_strain;					\
      stress_val += size_strain;					\
    }									\
  }                                                                     \


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent)		\
  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(el_type);	\
  UInt size_strain          = spatial_dimension * spatial_dimension;	\
  									\
  UInt nb_element;							\
  Real * strain_val;							\
  Real * tangent_val;							\
  									\
  if(ghost_type == _not_ghost) {					\
    nb_element   = element_filter[el_type]->getSize();			\
    stress[el_type]->resize(nb_element*nb_quadrature_points);		\
    strain_val = strain[el_type]->values;				\
  } else {								\
    nb_element = ghost_element_filter[el_type]->getSize();		\
    ghost_stress[el_type]->resize(nb_element*nb_quadrature_points);	\
    strain_val = ghost_strain[el_type]->values;				\
  }									\
  tangent_val = tangent.values;						\
  size_tangent = getTangentStiffnessVoigtSize();			\
  size_tangent *= size_tangent;						\
  									\
  if (nb_element == 0) return;						\
  									\
  for (UInt el = 0; el < nb_element; ++el) {				\
    for (UInt q = 0; q < nb_quadrature_points; ++q) {			\


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END			\
      strain_val += size_strain;					\
      tangent_val += size_tangent;					\
    }									\
  }                                                                     \


#endif /* __AKANTU_MATERIAL_HH__ */
