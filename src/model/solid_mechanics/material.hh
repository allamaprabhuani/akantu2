/**
 * @file   material.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
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

  /// assemble the residual for this material
  void assembleResidual(GhostType ghost_type);

  /// compute the stresses for this material
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);
  virtual void computeAllNonLocalStresses(__attribute__((unused)) GhostType ghost_type = _not_ghost) {};
  virtual void computeAllStressesFromTangentModuli(GhostType ghost_type = _not_ghost);

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


  template<UInt dim>
  void computeAllStressesFromTangentModuli(const ElementType & type, GhostType ghost_type);

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
  inline void transferBMatrixToSymVoigtBMatrix(const types::RMatrix & B,
					       types::RMatrix & Bvoigt,
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
  inline void buildElementalFieldInterpolationCoodinates(const types::RMatrix & coordinates,
							 types::RMatrix & coordMatrix);

  /// get the size of the coordiante matrix used in the interpolation
  template <ElementType type>
  inline UInt getSizeElementalFieldInterpolationCoodinates();


  /* ------------------------------------------------------------------------ */
  /* Function for all materials                                               */
  /* ------------------------------------------------------------------------ */

protected:
  /// compute the potential energy for a quadrature point
  inline void computePotentialEnergyOnQuad(types::RMatrix & grad_u,
					   types::RMatrix & sigma,
					   Real & epot);

  /// compute the potential energy for an element
  virtual void computePotentialEnergyByElement(ElementType type, UInt index,
					       types::RVector & epot_on_quad_points);

protected:
  /// allocate an internal vector
  template<typename T>
  void initInternalVector(ByElementTypeVector<T> & vect,
			  UInt nb_component,
			  bool temporary = false,
			  ElementKind element_kind = _ek_regular);

public:
  /// resize an internal vector
  template<typename T>
  void resizeInternalVector(ByElementTypeVector<T> & vect,
			    ElementKind element_kind = _ek_regular) const;

  /* ------------------------------------------------------------------------ */
  template<UInt dim>
  inline void gradUToF(const types::RMatrix & grad_u, types::RMatrix & F);
  inline void rightCauchy(const types::RMatrix & F, types::RMatrix & C);
  inline void leftCauchy (const types::RMatrix & F, types::RMatrix & B);
  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  virtual inline UInt getNbDataForElements(const Vector<Element> & elements,
					   SynchronizationTag tag) const;

  virtual inline void packElementData(CommunicationBuffer & buffer,
				      const Vector<Element> & elements,
				      SynchronizationTag tag) const;

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
					const Vector<Element> & elements,
					SynchronizationTag tag);

  template<typename T>
  inline void packElementDataHelper(const ByElementTypeVector<T> & data_to_pack,
				    CommunicationBuffer & buffer,
				    const Vector<Element> & elements) const;

  template<typename T>
  inline void unpackElementDataHelper(ByElementTypeVector<T> & data_to_unpack,
				      CommunicationBuffer & buffer,
				      const Vector<Element> & elements) const;

protected:
  template<typename T, bool pack_helper>
  inline void packUnpackElementDataHelper(ByElementTypeVector<T> & data_to_pack,
					  CommunicationBuffer & buffer,
					  const Vector<Element> & elements) const;
public:
  inline UInt getNbQuadraturePoints(const Vector<Element> & elements) const;

public:
  /* ------------------------------------------------------------------------ */
  virtual inline void onElementsAdded(const Vector<Element> & element_list);

  virtual inline void onElementsRemoved(const Vector<Element> & element_list,
					const ByElementTypeUInt & new_numbering);

protected:
  template<typename T>
  void removeQuadraturePointsFromVectors(ByElementTypeVector<T> & data,
					 const ByElementTypeUInt & new_numbering);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Model, *model, const SolidMechanicsModel &)

  AKANTU_GET_MACRO(ID, id, const ID &);
  AKANTU_GET_MACRO(Rho, rho, Real);
  AKANTU_SET_MACRO(Rho, rho, Real);

  /// return the potential energy for the subset of elements contained by the material
  Real getPotentialEnergy();
  /// return the potential energy for the provided element
  Real getPotentialEnergy(ElementType & type, UInt index);

  /// return the energy (identified by id) for the subset of elements contained by the material
  virtual Real getEnergy(std::string energy_id);
  /// return the energy (identified by id) for the provided element
  virtual Real getEnergy(std::string energy_id, ElementType type, UInt index);


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
  Vector<Real>::iterator<types::RMatrix> strain_it =			\
    this->strain(el_type, ghost_type).begin(spatial_dimension,		\
					    spatial_dimension);		\
  Vector<Real>::iterator<types::RMatrix> strain_end =			\
    this->strain(el_type, ghost_type).end(spatial_dimension,		\
					  spatial_dimension);		\
  Vector<Real>::iterator<types::RMatrix> stress_it =			\
    this->stress(el_type, ghost_type).begin(spatial_dimension,		\
					    spatial_dimension);		\
  									\
  for(;strain_it != strain_end; ++strain_it, ++stress_it) {		\
    types::RMatrix & __attribute__((unused)) grad_u = *strain_it;	\
    types::RMatrix & __attribute__((unused)) sigma  = *stress_it

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END			\
  }									\


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_mat)	\
  Vector<Real>::iterator<types::RMatrix> strain_it =			\
    this->strain(el_type, ghost_type).begin(spatial_dimension,		\
					    spatial_dimension);		\
  Vector<Real>::iterator<types::RMatrix> strain_end =			\
    this->strain(el_type, ghost_type).end(spatial_dimension,		\
					  spatial_dimension);		\
  									\
  UInt tangent_size =							\
    this->getTangentStiffnessVoigtSize(spatial_dimension);		\
  Vector<Real>::iterator<types::RMatrix> tangent_it =			\
    tangent_mat.begin(tangent_size,					\
		      tangent_size);					\
  									\
  for(;strain_it != strain_end; ++strain_it, ++tangent_it) {		\
    types::RMatrix & __attribute__((unused)) grad_u  = *strain_it;	\
    types::RMatrix & tangent = *tangent_it


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END			\
  }


/* -------------------------------------------------------------------------- */
#define INSTANSIATE_MATERIAL(mat_name)			\
  template class mat_name<1>;				\
  template class mat_name<2>;				\
  template class mat_name<3>

/* -------------------------------------------------------------------------- */
/* Material list                                                              */
/* -------------------------------------------------------------------------- */
// elastic materials
#include "material_elastic.hh"

#define AKANTU_CORE_MATERIAL_LIST					\
  ((2, (elastic            , MaterialElastic                      )))


#if defined(AKANTU_EXTRA_MATERIALS)
#  include "material_extra_includes.hh"
#else
#  define AKANTU_EXTRA_MATERIAL_LIST
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
#  include "material_cohesive_includes.hh"
#else
#  define AKANTU_COHESIVE_MATERIAL_LIST
#endif

#if defined(AKANTU_DAMAGE_NON_LOCAL)
#  include "material_non_local_includes.hh"
#else
#  define AKANTU_DAMAGE_NON_LOCAL_MATERIAL_LIST
#endif

#define AKANTU_MATERIAL_LIST			\
  AKANTU_CORE_MATERIAL_LIST			\
  AKANTU_EXTRA_MATERIAL_LIST			\
  AKANTU_COHESIVE_MATERIAL_LIST			\
  AKANTU_DAMAGE_NON_LOCAL_MATERIAL_LIST


#endif /* __AKANTU_MATERIAL_HH__ */
