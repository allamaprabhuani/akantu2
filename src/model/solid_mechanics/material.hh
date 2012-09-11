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
#include "parser.hh"
#include "data_accessor.hh"
#include "material_parameters.hh"

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
class Material : protected Memory, public DataAccessor, public Parsable, public MeshEventHandler {
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

  /// register a material parameter
  template<class T>
  void registerParam(std::string name, T & variable, T default_value,
		     ParamAccessType type,
		     std::string description = "");

  template<class T>
  void registerParam(std::string name, T & variable, ParamAccessType type,
		     std::string description = "");


  template<typename T>
  void registerInternal(__attribute__((unused)) ByElementTypeVector<T> & vect) { AKANTU_DEBUG_TO_IMPLEMENT(); }

  /// read parameter from file
  virtual bool parseParam(const std::string & key, const std::string & value,
			  const ID & id);

  /// function called to update the internal parameters when the modifiable
  /// parameters are modified
  virtual void updateInternalParameters() {}

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// compute the residual for this material
  virtual void updateResidual(GhostType ghost_type = _not_ghost);


  void assembleResidual(GhostType ghost_type);

  /// compute the residual for this material
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);
  virtual void computeAllNonLocalStresses(__attribute__((unused)) GhostType ghost_type = _not_ghost) {};

  /// set material to steady state
  void setToSteadyState(GhostType ghost_type = _not_ghost);

  /// compute the stiffness matrix
  virtual void assembleStiffnessMatrix(GhostType ghost_type);

  /// compute the stable time step for an element of size h
  virtual Real getStableTimeStep(__attribute__((unused)) Real h,
				 __attribute__((unused)) const Element & element = ElementNull)  {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed() const { AKANTU_DEBUG_TO_IMPLEMENT(); };

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed() const { AKANTU_DEBUG_TO_IMPLEMENT(); };

  /// add an element to the local mesh filter
  inline UInt addElement(const ElementType & type,
			 UInt element,
			 const GhostType & ghost_type);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /**
   * interpolate stress on given positions for each element by means
   * of a geometrical interpolation on quadrature points
   */
  virtual void interpolateStress(const ElementType type,
				 Vector<Real> & result);

  /**
   * function to initialize the elemental field interpolation
   * function by inverting the quadrature points' coordinates
   */
  virtual void initElementalFieldInterpolation(ByElementTypeReal & interpolation_points_coordinates);

protected:

  /// constitutive law
  virtual void computeStress(__attribute__((unused)) ElementType el_type,
			     __attribute__((unused)) GhostType ghost_type = _not_ghost)  {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }


  /// set the material to steady state (to be implemented for materials that need it)
  virtual void setToSteadyState(__attribute__((unused)) ElementType el_type,
				__attribute__((unused)) GhostType ghost_type = _not_ghost) {}

  /// compute the tangent stiffness matrix
  virtual void computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
				    __attribute__((unused)) Vector<Real> & tangent_matrix,
				    __attribute__((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute the potential energy
  virtual void computePotentialEnergy(ElementType el_type,
				      GhostType ghost_type = _not_ghost);

  template<UInt dim>
  void assembleStiffnessMatrix(const ElementType & type,
			       GhostType ghost_type);

  /// transfer the B matrix to a Voigt notation B matrix
  template<UInt dim>
  inline void transferBMatrixToSymVoigtBMatrix(const types::Matrix & B,
					       types::Matrix & Bvoigt,
					       UInt nb_nodes_per_element) const;

  inline UInt getTangentStiffnessVoigtSize(UInt spatial_dimension) const;


  /// compute the potential energy by element
  void computePotentialEnergyByElement();

  /// compute the coordinates of the quadrature points
  void computeQuadraturePointsCoordinates(ByElementTypeReal & quadrature_points_coordinates,
					  const GhostType & ghost_type) const;

  /// interpolate an elemental field on given points for each element
  template <ElementType type>
  void interpolateElementalField(const Vector<Real> & field,
				 Vector<Real> & result);

  /// template function to initialize the elemental field interpolation
  template <ElementType type>
  void initElementalFieldInterpolation(const Vector<Real> & quad_coordinates,
				       const Vector<Real> & interpolation_points_coordinates);

  /// build the coordinate matrix for the interpolation on elemental field
  template <ElementType type>
  inline void buildElementalFieldInterpolationCoodinates(const types::Matrix & coordinates,
							 types::Matrix & coordMatrix);

  /// get the size of the coordiante matrix used in the interpolation
  template <ElementType type>
  inline UInt getSizeElementalFieldInterpolationCoodinates();

  /// extract the field values corresponding to the quadrature points used for the interpolation
  template <ElementType type>
  inline void extractElementalFieldForInterplation(const Vector<Real> & field,
						   Vector<Real> & filtered_field);

  /* ------------------------------------------------------------------------ */
  /* Function for all materials                                               */
  /* ------------------------------------------------------------------------ */
protected:
  /// compute the potential energy for on element
  inline void computePotentialEnergyOnQuad(types::Matrix & grad_u,
					   types::Matrix & sigma,
					   Real & epot);

protected:
  /// allocate an internal vector
  template<typename T>
  void initInternalVector(ByElementTypeVector<T> & vect,
			  UInt nb_component,
			  ElementKind element_kind = _ek_regular);

public:
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

  virtual inline UInt getNbDataToPack  (const Element & element, SynchronizationTag tag) const;
  virtual inline UInt getNbDataToUnpack(const Element & element, SynchronizationTag tag) const;

  virtual UInt getNbDataToPack(__attribute__((unused)) SynchronizationTag tag) const {
    return 0;
  }

  virtual UInt getNbDataToUnpack(__attribute__((unused)) SynchronizationTag tag) const {
    return 0;
  }


  virtual inline void packData(CommunicationBuffer & buffer,
			       const Element & element,
			       SynchronizationTag tag) const;

  virtual void packData(__attribute__((unused)) CommunicationBuffer & buffer,
			__attribute__((unused)) const UInt index,
			__attribute__((unused)) SynchronizationTag tag) const {
  }

  virtual inline void unpackData(CommunicationBuffer & buffer,
				 const Element & element,
				 SynchronizationTag tag);

  virtual void unpackData(__attribute__((unused)) CommunicationBuffer & buffer,
			  __attribute__((unused)) const UInt index,
			  __attribute__((unused)) SynchronizationTag tag) {
  }

  /* ------------------------------------------------------------------------ */
  virtual inline void onElementsAdded(const Vector<Element> & element_list);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Model, *model, const SolidMechanicsModel &)

  AKANTU_GET_MACRO(ID, id, const ID &);
  AKANTU_GET_MACRO(Rho, rho, Real);
  AKANTU_SET_MACRO(Rho, rho, Real);

  Real getPotentialEnergy();
  virtual Real getEnergy(std::string type);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Strain, strain, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PotentialEnergy, potential_energy, Real);

  bool isNonLocal() const { return is_non_local; }

  const Vector<Real> & getVector(const ID & id, const ElementType & type, const GhostType & ghost_type = _not_ghost) const;
  Vector<Real> & getVector(const ID & id, const ElementType & type, const GhostType & ghost_type = _not_ghost);

  template<typename T>
  inline T getParam(const ID & param) const;

  template<typename T>
  inline void setParam(const ID & param, T value);

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

  /// elemental field interpolation coordinates
  ByElementTypeReal interpolation_inverse_coordinates;

  /// elemental field interpolation points
  ByElementTypeReal interpolation_points_matrices;

  /// list of the paramters
  MaterialParameters params;

private:
  /// boolean to know if the material has been initialized
  bool is_init;

  std::map<ID, ByElementTypeReal *> internal_vectors_real;
  std::map<ID, ByElementTypeUInt *> internal_vectors_uint;
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

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type) \
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
    types::Matrix & __attribute__((unused)) grad_u  = *strain_it;	\
    types::Matrix & tangent = *tangent_it


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END			\
  }


/* -------------------------------------------------------------------------- */
/* Material list                                                              */
/* -------------------------------------------------------------------------- */
#define AKANTU_CORE_MATERIAL_LIST					\
  ((2, (elastic            , MaterialElastic                      )))	\
  ((2, (elastic_orthotropic, MaterialElasticOrthotropic           )))	\
  ((2, (neohookean         , MaterialNeohookean                   )))	\
  ((2, (sls_deviatoric     , MaterialStandardLinearSolidDeviatoric)))	\
  ((2, (ve_stiffness_prop  , MaterialStiffnessProportional        )))


#define AKANTU_COHESIVE_MATERIAL_LIST					\
  ((2, (cohesive_bilinear      , MaterialCohesiveBilinear     )))	\
  ((2, (cohesive_linear        , MaterialCohesiveLinear       )))	\
  ((2, (cohesive_linear_extrinsic, MaterialCohesiveLinearExtrinsic )))	\
  ((2, (cohesive_linear_exponential_extrinsic, MaterialCohesiveLinearExponentialExtrinsic ))) \
  ((2, (cohesive_exponential   , MaterialCohesiveExponential  )))


#define  AKANTU_DAMAGE_MATERIAL_LIST					\
  ((2, (damage_linear          , MaterialDamageLinear         )))	\
  ((2, (marigo                 , MaterialMarigo               )))	\
  ((2, (mazars                 , MaterialMazars               )))	\
  ((2, (vreepeerlings          , MaterialVreePeerlings        )))

#ifdef AKANTU_DAMAGE_NON_LOCAL
#define AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST			\
  ((stress_wf, StressBasedWeightFunction  ))				\
  ((damage_wf, DamagedWeightFunction      ))				\
  ((remove_wf, RemoveDamagedWeightFunction))				\
  ((base_wf,   BaseWeightFunction         ))


#define AKANTU_MATERIAL_VREEPEERLINGS_WEIGHT_FUNCTION_TMPL_LIST		\
  AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST				\
  ((removed_damrate_wf, RemoveDamagedWithDamageRateWeightFunction))

#define AKANTU_DAMAGE_NON_LOCAL_MATERIAL_LIST				\
  ((3, (marigo_non_local       , MaterialMarigoNonLocal,		\
	AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST)))			\
  ((2, (mazars_non_local       , MaterialMazarsNonLocal       )))	\
  ((3, (vreepeerlings_non_local, MaterialVreePeerlingsNonLocal,		\
	AKANTU_MATERIAL_VREEPEERLINGS_WEIGHT_FUNCTION_TMPL_LIST)))

#else
#  define AKANTU_DAMAGE_NON_LOCAL_MATERIAL_LIST
#endif

#define AKANTU_MATERIAL_LIST			\
  AKANTU_CORE_MATERIAL_LIST			\
  AKANTU_COHESIVE_MATERIAL_LIST			\
  AKANTU_DAMAGE_MATERIAL_LIST			\
  AKANTU_DAMAGE_NON_LOCAL_MATERIAL_LIST


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
#include "material_neohookean.hh"
#include "material_elastic_orthotropic.hh"

// visco-elastic materials
#include "material_stiffness_proportional.hh"
#include "material_standard_linear_solid_deviatoric.hh"


// damage materials
#include "material_damage.hh"
#include "material_marigo.hh"
#include "material_mazars.hh"
#include "material_damage_linear.hh"
#include "material_vreepeerlings.hh"

#if defined(AKANTU_DAMAGE_NON_LOCAL)
#  include "material_marigo_non_local.hh"
#  include "material_mazars_non_local.hh"
#  include "material_vreepeerlings_non_local.hh"
#endif

// cohesive materials
#include "material_cohesive.hh"
#include "material_cohesive_linear.hh"
#include "material_cohesive_bilinear.hh"
#include "material_cohesive_linear_extrinsic.hh"
#include "material_cohesive_exponential.hh"
#include "material_cohesive_linear_exponential_extrinsic.hh"


#if defined(__INTEL_COMPILER)
#pragma warning ( pop )
#endif //defined(__INTEL_COMPILER)


#endif /* __AKANTU_MATERIAL_HH__ */
