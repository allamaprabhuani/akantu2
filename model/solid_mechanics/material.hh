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
  class Model;
  class SolidMechanicsModel;
  class CommunicationBuffer;
}

__BEGIN_AKANTU__

/** 
 * Interface of all materials
 * Prerequisites for a new material
 * - inherit from this class
 * - implement the following methods:
 * \code
 *  virtual Real getStableTimeStep(Real h, const Element & element = ElementNull);
 *
 *  virtual void computeStress(ElementType el_type,
 *                             GhostType ghost_type = _not_ghost);
 *
 *  virtual void computeTangentStiffness(const ElementType & el_type,
 *                                       Vector<Real> & tangent_matrix,
 *                                       GhostType ghost_type = _not_ghost);
 * \endcode
 *
 */
class Material : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Material(Model & model, const ID & id = "");
  virtual ~Material();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// read properties
  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// compute the residual for this material
  void updateResidual(Vector<Real> & current_position,
		      GhostType ghost_type = _not_ghost);


  /// compute the stiffness matrix
  void assembleStiffnessMatrix(Vector<Real> & current_position,
			       GhostType ghost_type);

  /// compute the stable time step for an element of size h
  virtual Real getStableTimeStep(Real h, const Element & element = ElementNull) = 0;

  /// add an element to the local mesh filter
  inline void addElement(const ElementType & type,
			 UInt element,
			 const GhostType & ghost_type);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

protected:

  /// constitutive law
  virtual void computeStress(ElementType el_type,
			     GhostType ghost_type = _not_ghost) = 0;

  /// constitutive law
  virtual void computeNonLocalStress(ElementType el_type,
				     GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  /// compute the tangent stiffness matrix
  virtual void computeTangentStiffness(const ElementType & el_type,
				       Vector<Real> & tangent_matrix,
				       GhostType ghost_type = _not_ghost) = 0;

  /// compute the potential energy
  virtual void computePotentialEnergy(ElementType el_type,
				      GhostType ghost_type = _not_ghost);

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
  /// compute the potential energy for on element
  inline void computePotentialEnergy(Real * F, Real * sigma, Real * epot);

  /// allocate an internal vector
  template<typename T>
  void initInternalVector(ByElementTypeVector<T> & vect,
			  UInt nb_component);

  /// resize an internal vector
  template<typename T>
  void resizeInternalVector(ByElementTypeVector<T> & vect);

  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
public:

  inline virtual UInt getNbDataToPack(const Element & element,
				      SynchronizationTag tag);

  inline virtual UInt getNbDataToUnpack(const Element & element,
					SynchronizationTag tag);

  inline virtual void packData(CommunicationBuffer & buffer,
			       const Element & element,
			       SynchronizationTag tag);

  inline virtual void unpackData(CommunicationBuffer & buffer,
				 const Element & element,
				 SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(ID, id, const ID &);
  AKANTU_GET_MACRO(Rho, rho, Real);

  Real getPotentialEnergy();

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Strain, strain, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PotentialEnergy, potential_energy, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the material
  ID id;

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

  /// is the vector for potential energy initialized
  //  bool potential_energy_vector;

  /// potential energy by element
  ByElementTypeReal potential_energy;

  /// boolean to know if the material has been initialized
  bool is_init;

  /// tell if using in non local mode or not
  bool is_non_local;

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

/* -------------------------------------------------------------------------- */
/* Auto loop                                                                  */
/* -------------------------------------------------------------------------- */

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN			\
  UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(el_type, ghost_type); \
  UInt size_strain = spatial_dimension * spatial_dimension;		\
  									\
  UInt nb_element = element_filter(el_type, ghost_type).getSize();	\
  if (nb_element == 0) return;						\
									\
  Vector<Real> & stress_tmp = stress(el_type, ghost_type);		\
  stress_tmp.resize(nb_element*nb_quadrature_points);			\
  Vector<Real> & strain_tmp = strain(el_type, ghost_type);		\
  									\
  Real * strain_val = strain_tmp.storage();				\
  Real * stress_val = stress_tmp.storage();				\
  									\
  for (UInt el = 0; el < nb_element; ++el) {				\
    for (UInt q = 0; q < nb_quadrature_points; ++q) {			\


#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END			\
      strain_val += size_strain;					\
      stress_val += size_strain;					\
    }									\
  }									\


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent)		\
  UInt nb_quadrature_points =						\
    model->getFEM().getNbQuadraturePoints(el_type, ghost_type);		\
  UInt size_strain = spatial_dimension * spatial_dimension;		\
									\
  UInt nb_element = element_filter(el_type, ghost_type).getSize();	\
  if (nb_element == 0) return;						\
									\
  Vector<Real> & strain_tmp = strain(el_type, ghost_type);		\
  									\
  Real * strain_val = strain_tmp.storage();				\
									\
  Real * tangent_val = tangent.values;					\
  UInt size_tangent = getTangentStiffnessVoigtSize(spatial_dimension);	\
  size_tangent *= size_tangent;						\
  									\
  									\
  for (UInt el = 0; el < nb_element; ++el) {				\
    for (UInt q = 0; q < nb_quadrature_points; ++q) {			\


// Vector<Real> * stress_tmp = stress(el_type, ghost_type);
// stress_tmp->resize(nb_element*nb_quadrature_points);
// Real * stress_val = stress_tmp->values;


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END			\
      strain_val += size_strain;					\
      tangent_val += size_tangent;					\
    }									\
  }                                                                     \


/* -------------------------------------------------------------------------- */
/* Material list                                                              */
/* -------------------------------------------------------------------------- */

#define AKANTU_MATERIAL_LIST						\
  ((elastic        , MaterialElastic       ))				\
  ((elastic_caughey, MaterialElasticCaughey))				\
  ((damage         , MaterialDamage        ))				\
  ((mazars         , MaterialMazars        ))				\
  ((neohookean     , MaterialNeohookean    ))				\
  ((non_local      , MaterialNonLocal      ))

#include "materials/material_elastic.hh"
#include "materials/material_elastic_caughey.hh"
#include "materials/material_damage.hh"
#include "materials/material_mazars.hh"
#include "materials/material_neohookean.hh"
#include "materials/material_non_local.hh"

#endif /* __AKANTU_MATERIAL_HH__ */
