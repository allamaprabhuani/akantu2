/**
 * @file   material.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Jul 27 18:15:37 2010
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
#include "parsable.hh"
#include "data_accessor.hh"
#include "internal_field.hh"
#include "random_internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_HH__
#define __AKANTU_MATERIAL_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
  class Model;
  class SolidMechanicsModel;
}

__BEGIN_AKANTU__

template <UInt dim>
class VoigtHelper {
public:
  const static UInt size;
  // matrix of vector index I as function of tensor indices i,j
  const static UInt mat[dim][dim];
  // array of matrix indices ij as function of vector index I
  const static UInt vec[dim*dim][2];
  // factors to multiply the strain by for voigt notation
  const static Real factors[dim*(dim-(dim-1)/2)];
};
template <UInt dim> const UInt VoigtHelper<dim>::size = dim*(dim-(dim-1)/2);


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
 *                                       Array<Real> & tangent_matrix,
 *                                       GhostType ghost_type = _not_ghost);
 * \endcode
 *
 */
class Material : public Memory, public DataAccessor, public Parsable, public MeshEventHandler {
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

  template<typename T>
  void registerInternal(__attribute__((unused)) InternalField<T> & vect) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  template<typename T>
  void unregisterInternal(__attribute__((unused)) InternalField<T> & vect) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// function called to update the internal parameters when the modifiable
  /// parameters are modified
  virtual void updateInternalParameters() {}

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// compute the residual for this material
  virtual void updateResidual(GhostType ghost_type = _not_ghost);

  /// assemble the residual for this material
  virtual void assembleResidual(GhostType ghost_type);

  /// save the stress in the precious_stress if needed
  virtual void savePreviousState(GhostType ghost_type);

  /// compute the stresses for this material
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);
  virtual void computeAllNonLocalStresses(__attribute__((unused)) GhostType ghost_type = _not_ghost) {};
  virtual void computeAllStressesFromTangentModuli(GhostType ghost_type = _not_ghost);
  virtual void computeAllCauchyStresses(GhostType ghost_type = _not_ghost);

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
  virtual void interpolateStress(ByElementTypeReal & result,
				 const GhostType ghost_type = _not_ghost);

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

  /// 
  template<UInt dim>
  void assembleResidual(GhostType ghost_type);

  /// Computation of Cauchy stress tensor in the case of finite deformation
  template<UInt dim>
  void computeCauchyStress(__attribute__((unused)) ElementType el_type,
                          __attribute__((unused)) GhostType ghost_type = _not_ghost);

  template<UInt dim >
  inline void computeCauchyStressOnQuad(const Matrix<Real> & F, const Matrix<Real> & S,
        Matrix<Real> & cauchy);

  template<UInt dim>
  void computeAllStressesFromTangentModuli(const ElementType & type, GhostType ghost_type);

  /// set the material to steady state (to be implemented for materials that need it)
  virtual void setToSteadyState(__attribute__((unused)) ElementType el_type,
                                __attribute__((unused)) GhostType ghost_type = _not_ghost) {}

  /// compute the tangent stiffness matrix
  virtual void computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
                                    __attribute__((unused)) Array<Real> & tangent_matrix,
                                    __attribute__((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute the potential energy
  virtual void computePotentialEnergy(ElementType el_type,
                                      GhostType ghost_type = _not_ghost);

  template<UInt dim>
  void assembleStiffnessMatrix(const ElementType & type,
                               GhostType ghost_type);

  /// assembling in finite deformation
  template<UInt dim>
  void assembleStiffnessMatrixNL(const ElementType & type,
                                 GhostType ghost_type);

  template<UInt dim>
  void assembleStiffnessMatrixL2(const ElementType & type,
                                 GhostType ghost_type);


  /// transfer the B matrix to a Voigt notation B matrix
  template<UInt dim>
  inline void transferBMatrixToSymVoigtBMatrix(const Matrix<Real> & B,
                                               Matrix<Real> & Bvoigt,
                                               UInt nb_nodes_per_element) const;

  /// transfer the BNL and BL2 matrix to a Voigt notation B matrix (See Bathe et al. IJNME vol 9, 1975)
  template<UInt dim>
  inline void transferBMatrixToBNL(const Matrix<Real> & B,
                                   Matrix<Real> & Bvoigt,
                                   UInt nb_nodes_per_element) const;

  template<UInt dim>
  inline void transferBMatrixToBL2(const Matrix<Real> & B, const Matrix<Real> & grad_u,
                                   Matrix<Real> & Bvoigt,
                                   UInt nb_nodes_per_element) const;

  /// write the stress tensor in the Voigt notation.
  template<UInt dim>
  inline void SetCauchyStressArray(const Matrix<Real> & S_t, Matrix<Real> & Stress_vect);


  inline UInt getTangentStiffnessVoigtSize(UInt spatial_dimension) const;

  /// Size of the Stress matrix for the case of finite deformation see: Bathe et al, IJNME, Vol 9, 353-386, 1975
  inline UInt getCauchyStressMatrixSize(UInt spatial_dimension) const;

  /// Sets the stress matrix according to Bathe et al, IJNME, Vol 9, 353-386, 1975
  template<UInt dim>
  inline void SetCauchyStressMatrix(const Matrix<Real> & S_t, Matrix<Real> & Stress_matrix);

  /// compute the potential energy by element
  void computePotentialEnergyByElement();

public:
  /// compute the coordinates of the quadrature points
  void computeQuadraturePointsCoordinates(ByElementTypeReal & quadrature_points_coordinates,
                                          const GhostType & ghost_type) const;

protected:
  /// interpolate an elemental field on given points for each element
  template <ElementType type>
  void interpolateElementalField(const Array<Real> & field,
				 Array<Real> & result,
				 const GhostType ghost_type);

  /// template function to initialize the elemental field interpolation
  template <ElementType type>
  void initElementalFieldInterpolation(const Array<Real> & quad_coordinates,
				       const Array<Real> & interpolation_points_coordinates,
				       const GhostType ghost_type);

  /// build the coordinate matrix for the interpolation on elemental field
  template <ElementType type>
  inline void buildElementalFieldInterpolationCoodinates(const Matrix<Real> & coordinates,
                                                         Matrix<Real> & coordMatrix);

  /// build interpolation coordinates for basic linear elements
  inline void buildElementalFieldInterpolationCoodinatesLinear(const Matrix<Real> & coordinates,
							       Matrix<Real> & coordMatrix);

  /// build interpolation coordinates for basic quadratic elements
  inline void buildElementalFieldInterpolationCoodinatesQuadratic(const Matrix<Real> & coordinates,
								  Matrix<Real> & coordMatrix);

  /// get the size of the coordiante matrix used in the interpolation
  template <ElementType type>
  inline UInt getSizeElementalFieldInterpolationCoodinates(GhostType ghost_type = _not_ghost);


  /* ------------------------------------------------------------------------ */
  /* Function for all materials                                               */
  /* ------------------------------------------------------------------------ */

protected:
  /// compute the potential energy for a quadrature point
  virtual inline void computePotentialEnergyOnQuad(Matrix<Real> & grad_u,
                                           Matrix<Real> & sigma,
                                           Real & epot);

  /// compute the potential energy for an element
  virtual void computePotentialEnergyByElement(ElementType type, UInt index,
                                               Vector<Real> & epot_on_quad_points);


public:
  /* ------------------------------------------------------------------------ */
  template<UInt dim>
  inline void gradUToF(const Matrix<Real> & grad_u, Matrix<Real> & F);
  inline void rightCauchy(const Matrix<Real> & F, Matrix<Real> & C);
  inline void leftCauchy (const Matrix<Real> & F, Matrix<Real> & B);

  template<UInt dim>
  inline void deformationJacobian(const Matrix<Real> & F, Real & J);

  inline void gradUToEpsilon(const Matrix<Real> & grad_u, Matrix<Real> & epsilon);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
                                           SynchronizationTag tag) const;

  virtual inline void packElementData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      SynchronizationTag tag) const;

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        SynchronizationTag tag);

  template<typename T>
  inline void packElementDataHelper(const ByElementTypeArray<T> & data_to_pack,
                                    CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const ID & fem_id = ID()) const;

  template<typename T>
  inline void unpackElementDataHelper(ByElementTypeArray<T> & data_to_unpack,
                                      CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      const ID & fem_id = ID());

public:
  /* ------------------------------------------------------------------------ */
  virtual inline void onElementsAdded(const Array<Element> & element_list,
                                      const NewElementsEvent & event);

  virtual inline void onElementsRemoved(const Array<Element> & element_list,
                                        const ByElementTypeUInt & new_numbering,
                                        const RemovedElementsEvent & event);

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
  AKANTU_GET_MACRO(Strain, strain, const ByElementTypeReal &);
  AKANTU_GET_MACRO(Stress, stress, const ByElementTypeReal &);
  AKANTU_GET_MACRO(ElementFilter, element_filter, const ByElementTypeUInt &);

  bool isNonLocal() const { return is_non_local; }

  const Array<Real> & getArray(const ID & id, const ElementType & type, const GhostType & ghost_type = _not_ghost) const;
  Array<Real> & getArray(const ID & id, const ElementType & type, const GhostType & ghost_type = _not_ghost);

  const ByElementTypeArray<Real> & getInternal(const ID & id) const;
  ByElementTypeArray<Real> & getInternal(const ID & id);

  bool isFiniteDeformation() const { return finite_deformation; }
  bool isInelasticDeformation() const { return inelastic_deformation; }

protected:

  bool isInit() const { return is_init; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// boolean to know if the material has been initialized
  bool is_init;

  std::map<ID, InternalField<Real> *> internal_vectors_real;
  std::map<ID, InternalField<UInt> *> internal_vectors_uint;

protected:
  /// id of the material
  ID id;

  /// Finite deformation
  bool finite_deformation;

  /// Finite deformation
  bool inelastic_deformation;


  /// material name
  std::string name;

  /// The model to witch the material belong
  SolidMechanicsModel * model;

  /// density : rho
  Real rho;

  /// spatial dimension
  UInt spatial_dimension;

  /// list of element handled by the material
  ByElementTypeArray<UInt> element_filter;

  /// stresses arrays ordered by element types
  InternalField<Real> stress;

  /// strains arrays ordered by element types
  InternalField<Real> strain;

  /// stress increment arrays ordered by element types (Finite deformation)
  InternalField<Real> delta_stress;

  /// strain increment arrays ordered by element types (Finite deformation)
  InternalField<Real> delta_strain;

  /// strain increment arrays ordered by element types (inelastic deformation)
  InternalField<Real> inelas_strain;

  /// Second Piola-Kirchhoff stress tensor arrays ordered by element types (Finite deformation)
  InternalField<Real> piola_kirchhoff_stress;

  /// potential energy by element
  InternalField<Real> potential_energy;

  /// tell if using in non local mode or not
  bool is_non_local;

  /// tell if the material need the previous stress state
  bool use_previous_stress;

  /// previous stresses
  InternalField<Real> previous_stress;

  /// tell if the material need the previous strain state
  bool use_previous_strain;

  /// previous strain
  InternalField<Real> previous_strain;

  /// elemental field interpolation coordinates
  InternalField<Real> interpolation_inverse_coordinates;

  /// elemental field interpolation points
  InternalField<Real> interpolation_points_matrices;
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

#include "internal_field_tmpl.hh"
#include "random_internal_field_tmpl.hh"

/* -------------------------------------------------------------------------- */
/* Auto loop                                                                  */
/* -------------------------------------------------------------------------- */

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type) \
  Array<Real>::iterator< Matrix<Real> > strain_it =			\
    this->strain(el_type, ghost_type).begin(spatial_dimension,		\
                                            spatial_dimension);		\
  Array<Real>::iterator< Matrix<Real> > strain_end =			\
    this->strain(el_type, ghost_type).end(spatial_dimension,		\
                                          spatial_dimension);		\
                                                                        \
  this->stress(el_type,                                                 \
               ghost_type).resize(this->strain(el_type,                 \
                                               ghost_type).getSize());  \
                                                                        \
  Array<Real>::iterator< Matrix<Real> > stress_it =			\
    this->stress(el_type, ghost_type).begin(spatial_dimension,		\
                                            spatial_dimension);		\
                                                                        \
  if(this->isFiniteDeformation())                                       \
    Array<Real>::iterator< Matrix<Real> > stress_it =			\
      this->piola_kirchhoff_stress(el_type, ghost_type).begin(spatial_dimension,        \
                                              spatial_dimension);       \
  for(;strain_it != strain_end; ++strain_it, ++stress_it) {		\
    Matrix<Real> & __attribute__((unused)) grad_u = *strain_it;         \
    Matrix<Real> & __attribute__((unused)) sigma  = *stress_it

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END			\
  }									\


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_mat)	\
  Array<Real>::iterator< Matrix<Real> > strain_it =			\
    this->strain(el_type, ghost_type).begin(spatial_dimension,		\
                                            spatial_dimension);		\
  Array<Real>::iterator< Matrix<Real> > strain_end =			\
    this->strain(el_type, ghost_type).end(spatial_dimension,		\
                                          spatial_dimension);		\
  Array<Real>::iterator< Matrix<Real> > sigma_it =                      \
    this->stress(el_type, ghost_type).begin(spatial_dimension,          \
                                                  spatial_dimension);   \
                                                                        \
  tangent_mat.resize(this->strain(el_type, ghost_type).getSize());      \
                                                                        \
  UInt tangent_size =							\
    this->getTangentStiffnessVoigtSize(spatial_dimension);		\
  Array<Real>::iterator< Matrix<Real> > tangent_it =			\
    tangent_mat.begin(tangent_size,					\
                      tangent_size);					\
                                                                        \
  for(;strain_it != strain_end; ++strain_it, ++sigma_it, ++tangent_it) { \
    Matrix<Real> & __attribute__((unused)) grad_u  = *strain_it;	\
    Matrix<Real> & __attribute__((unused)) sigma_tensor = *sigma_it;    \
    Matrix<Real> & tangent = *tangent_it


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END			\
  }                                                                     \


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_FINITE_DEFORMATION_BEGIN(tangent_mat)	\
  Array<Real>::iterator< Matrix<Real> > strain_it =			\
    this->strain(el_type, ghost_type).begin(spatial_dimension,		\
                                            spatial_dimension);		\
  Array<Real>::iterator< Matrix<Real> > strain_end =			\
    this->strain(el_type, ghost_type).end(spatial_dimension,		\
                                          spatial_dimension);		\
  Array<Real>::iterator< Matrix<Real> > green_it =                      \
    this->delta_strain(el_type, ghost_type).begin(spatial_dimension,    \
                                                  spatial_dimension);   \
  Array<Real>::iterator< Matrix<Real> > sigma_it =                      \
    this->stress(el_type, ghost_type).begin(spatial_dimension,          \
                                                  spatial_dimension);   \
                                                                        \
  tangent_mat.resize(this->strain(el_type, ghost_type).getSize());      \
                                                                        \
  UInt tangent_size =							\
    this->getTangentStiffnessVoigtSize(spatial_dimension);		\
  Array<Real>::iterator< Matrix<Real> > tangent_it =			\
    tangent_mat.begin(tangent_size,					\
                      tangent_size);					\
                                                                        \
  for(;strain_it != strain_end; ++strain_it, ++green_it, ++sigma_it, ++tangent_it) { \
    Matrix<Real> & __attribute__((unused)) grad_u  = *strain_it;	\
    Matrix<Real> & __attribute__((unused)) grad_delta_u = *green_it;    \
    Matrix<Real> & __attribute__((unused)) sigma_tensor = *sigma_it;    \
    Matrix<Real> & tangent = *tangent_it


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_FINITE_DEFORMATION_END	\
  }                                                                     \


/* -------------------------------------------------------------------------- */
#define INSTANSIATE_MATERIAL(mat_name)			\
  template class mat_name<1>;				\
  template class mat_name<2>;				\
  template class mat_name<3>

#endif /* __AKANTU_MATERIAL_HH__ */
