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
#include "aka_common.hh"
#include "aka_memory.hh"
//#include "fem.hh"
//#include "mesh.hh"
#include "data_accessor.hh"
//#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_HH__
#define __AKANTU_MATERIAL_HH__

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
class Material : protected Memory, public DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Material(SolidMechanicsModel & model, const ID & id = "");
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
  virtual void updateResidual(Vector<Real> & displacement,
			      GhostType ghost_type = _not_ghost);

  void assembleResidual(GhostType ghost_type);

  /// compute the residual for this material
  virtual void computeStress(Vector<Real> & current_position,
			     GhostType ghost_type = _not_ghost);

  /// set material to steady state
  void setToSteadyState(GhostType ghost_type = _not_ghost);

  /// compute the stiffness matrix
  void assembleStiffnessMatrix(Vector<Real> & current_position,
			       GhostType ghost_type);

  /// compute the stable time step for an element of size h
  virtual Real getStableTimeStep(Real h, const Element & element = ElementNull) = 0;

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed() { AKANTU_DEBUG_TO_IMPLEMENT(); };

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed() { AKANTU_DEBUG_TO_IMPLEMENT(); };

  /// add an element to the local mesh filter
  inline UInt addElement(const ElementType & type,
			 UInt element,
			 const GhostType & ghost_type);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:

  /// constitutive law
  virtual void computeStress(ElementType el_type,
			     GhostType ghost_type = _not_ghost) = 0;


  /// set the material to steady state (to be implemented for materials that need it)
  virtual void setToSteadyState(ElementType el_type,
				GhostType ghost_type = _not_ghost) {};

  // /// constitutive law
  // virtual void computeNonLocalStress(ElementType el_type,
  // 				     GhostType ghost_type = _not_ghost) {
  //   AKANTU_DEBUG_TO_IMPLEMENT();
  // };

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


  void computeQuadraturePointsCoordinates(ByElementTypeReal & quadrature_points_coordinates) const;

  /* ------------------------------------------------------------------------ */
  /* Function for all materials                                               */
  /* ------------------------------------------------------------------------ */
protected:
  /// compute the potential energy for on element
  inline void computePotentialEnergyOnQuad(types::Matrix & grad_u,
					   types::Matrix & sigma,
					   Real & epot);


public:
  /// allocate an internal vector
  template<typename T>
  void initInternalVector(ByElementTypeVector<T> & vect,
			  UInt nb_component,
			  ElementKind element_kind = _ek_regular) const;

  /// resize an internal vector
  template<typename T>
  void resizeInternalVector(ByElementTypeVector<T> & vect,
			    ElementKind element_kind = _ek_regular) const;

  /* ------------------------------------------------------------------------ */
  template<UInt dim>
  inline void gradUToF(const types::Matrix & grad_u, types::Matrix & F);
  inline void rightCauchy(const types::Matrix & F, types::Matrix & C);
  inline void leftCauchy (const types::Matrix & F, types::Matrix & B);
  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  virtual inline UInt getNbDataToPack(__attribute__((unused)) const Element & element,
				      __attribute__((unused)) SynchronizationTag tag) const {
    return 0;
  }

  virtual inline UInt getNbDataToUnpack(__attribute__((unused)) const Element & element,
					__attribute__((unused)) SynchronizationTag tag) const {
    return 0;
  }

  virtual UInt getNbDataToPack(__attribute__((unused)) SynchronizationTag tag) const {
    return 0;
  }

  virtual UInt getNbDataToUnpack(__attribute__((unused)) SynchronizationTag tag) const {
    return 0;
  }


  virtual inline void packData(__attribute__((unused)) CommunicationBuffer & buffer,
			       __attribute__((unused)) const Element & element,
			       __attribute__((unused)) SynchronizationTag tag) const {
  }

  virtual void packData(__attribute__((unused)) CommunicationBuffer & buffer,
			__attribute__((unused)) const UInt index,
                        __attribute__((unused)) SynchronizationTag tag) const {
  }

  virtual inline void unpackData(__attribute__((unused)) CommunicationBuffer & buffer,
				 __attribute__((unused)) const Element & element,
				 __attribute__((unused)) SynchronizationTag tag) {
  }

  virtual void unpackData(__attribute__((unused)) CommunicationBuffer & buffer,
			  __attribute__((unused)) const UInt index,
			  __attribute__((unused)) SynchronizationTag tag) {
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Model, *model, const SolidMechanicsModel &)

  AKANTU_GET_MACRO(ID, id, const ID &);
  AKANTU_GET_MACRO(Rho, rho, Real);
  AKANTU_SET_MACRO(Rho, rho, Real);

  Real getPotentialEnergy();

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Strain, strain, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PotentialEnergy, potential_energy, Real);

  const Vector<Real> & getVector(const ID & id, const ElementType & type, const GhostType & ghost_type = _not_ghost) const;

  virtual Real getParam(const ID & param) const;
  virtual void setParam(const ID & param, Real value);

protected:

  bool isInit() const { return is_init; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the material
  ID id;

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

  /// tell if using in non local mode or not
  bool is_non_local;

  /// spatial dimension
  UInt spatial_dimension;

private:
  /// boolean to know if the material has been initialized
  bool is_init;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "material_inline_impl.cc"
#endif

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
  Vector<Real>::iterator<types::Matrix> strain_it =			\
    this->strain(el_type, ghost_type).begin(spatial_dimension,		\
					    spatial_dimension);		\
  Vector<Real>::iterator<types::Matrix> strain_end =			\
    this->strain(el_type, ghost_type).end(spatial_dimension,		\
					  spatial_dimension);		\
  Vector<Real>::iterator<types::Matrix> stress_it =			\
    this->stress(el_type, ghost_type).begin(spatial_dimension,		\
					    spatial_dimension);		\
  									\
  for(;strain_it != strain_end; ++strain_it, ++stress_it) {		\
    types::Matrix & __attribute__((unused)) grad_u = *strain_it;	\
    types::Matrix & __attribute__((unused)) sigma  = *stress_it

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END			\
  }									\


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_mat)	\
  Vector<Real>::iterator<types::Matrix> strain_it =			\
    this->strain(el_type, ghost_type).begin(spatial_dimension,		\
					    spatial_dimension);		\
  Vector<Real>::iterator<types::Matrix> strain_end =			\
    this->strain(el_type, ghost_type).end(spatial_dimension,		\
					  spatial_dimension);		\
									\
  UInt tangent_size =							\
    this->getTangentStiffnessVoigtSize(spatial_dimension);		\
  Vector<Real>::iterator<types::Matrix> tangent_it =			\
    tangent_mat.begin(tangent_size,					\
		      tangent_size);					\
									\
  for(;strain_it != strain_end; ++strain_it, ++tangent_it) {		\
  types::Matrix & __attribute__((unused)) grad_u  = *strain_it;		\
    types::Matrix & tangent = *tangent_it


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END			\
  }


/* -------------------------------------------------------------------------- */
/* Material list                                                              */
/* -------------------------------------------------------------------------- */

#define AKANTU_MATERIAL_LIST						\
  ((elastic                , MaterialElastic              ))		\
  ((elastic_orthotropic    , MaterialElasticOrthotropic   ))		\
  ((viscoelastic           , MaterialViscoElastic         ))		\
  ((elastic_caughey        , MaterialElasticCaughey       ))		\
  ((neohookean             , MaterialNeohookean           ))		\
  ((damage_linear          , MaterialDamageLinear         ))		\
  ((marigo                 , MaterialMarigo               ))		\
  ((mazars                 , MaterialMazars               ))		\
  ((vreepeerlings          , MaterialVreePeerlings        ))		\
  ((marigo_non_local       , MaterialMarigoNonLocal       ))		\
  ((mazars_non_local       , MaterialMazarsNonLocal       ))		\
  ((vreepeerlings_non_local, MaterialVreePeerlingsNonLocal))		\
  ((cohesive_bilinear      , MaterialCohesiveBilinear     ))		\
  ((cohesive_linear        , MaterialCohesiveLinear       ))		\
  ((cohesive_linear_extrinsic, MaterialCohesiveLinearExtrinsic ))	\
  ((cohesive_linear_exponential_extrinsic, MaterialCohesiveLinearExponentialExtrinsic ))


#define INSTANSIATE_MATERIAL(mat_name)			\
  template class mat_name<1>;				\
  template class mat_name<2>;				\
  template class mat_name<3>

#if defined(__INTEL_COMPILER)
#pragma warning ( push )
/* warning #654: overloaded virtual function
   "akantu::Material::computeStress" is only partially overridden in
   class "akantu::Material*" */

#pragma warning ( disable : 654 )
#endif //defined(__INTEL_COMPILER)

/* -------------------------------------------------------------------------- */
// elastic materials
#include "material_elastic.hh"
#include "material_elastic_caughey.hh"
#include "material_viscoelastic.hh"
#include "material_neohookean.hh"
#include "material_elastic_orthotropic.hh"

// damage materials
#include "material_damage.hh"
#include "material_marigo.hh"
#include "material_mazars.hh"
#include "material_damage_linear.hh"
#include "material_vreepeerlings.hh"

#include "material_marigo_non_local.hh"
#include "material_mazars_non_local.hh"
#include "material_vreepeerlings_non_local.hh"

// cohesive materials
#include "material_cohesive.hh"
#include "material_cohesive_linear.hh"
#include "material_cohesive_bilinear.hh"
#include "material_cohesive_linear_extrinsic.hh"
#include "material_cohesive_linear_exponential_extrinsic.hh"


#if defined(__INTEL_COMPILER)
#pragma warning ( pop )
#endif //defined(__INTEL_COMPILER)


#endif /* __AKANTU_MATERIAL_HH__ */
