/**
 * @file   material.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Tue Sep 16 2014
 *
 * @brief  Mother class for all materials
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_voigthelper.hh"
#include "parser.hh"
#include "parsable.hh"
#include "data_accessor.hh"
#include "internal_field.hh"
#include "random_internal_field.hh"
#include "solid_mechanics_model_event_handler.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_HH__
#define __AKANTU_MATERIAL_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
  class Model;
  class SolidMechanicsModel;
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
 *                                       Array<Real> & tangent_matrix,
 *                                       GhostType ghost_type = _not_ghost);
 * \endcode
 *
 */

class Material : public Memory, public DataAccessor, public Parsable,
                 public MeshEventHandler,
                 protected SolidMechanicsModelEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
#if __cplusplus > 199711L
  Material(const Material & mat) = delete;
  Material & operator=(const Material & mat) = delete;
#endif

  /// Initialize material with defaults
  Material(SolidMechanicsModel & model, const ID & id = "");

  /// Initialize material with custom mesh & fe_engine
  Material(SolidMechanicsModel & model,
           UInt dim,
           const Mesh & mesh,
           FEEngine & fe_engine,
           const ID & id = "");

  /// Destructor
  virtual ~Material();

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Function that materials can/should reimplement                           */
  /* ------------------------------------------------------------------------ */
protected:
  /// constitutive law
  virtual void computeStress(__attribute__((unused)) ElementType el_type,
                             __attribute__((unused)) GhostType ghost_type = _not_ghost)  {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute the tangent stiffness matrix
  virtual void computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
                                    __attribute__((unused)) Array<Real> & tangent_matrix,
                                    __attribute__((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute the potential energy
  virtual void computePotentialEnergy(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the potential energy for an element
  virtual void computePotentialEnergyByElement(__attribute__((unused)) ElementType type,
                                               __attribute__((unused)) UInt index,
                                               __attribute__((unused)) Vector<Real> & epot_on_quad_points) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  virtual void updateEnergies(__attribute__((unused)) ElementType el_type,
                              __attribute__((unused)) GhostType ghost_type = _not_ghost) {  }

  virtual void updateEnergiesAfterDamage(__attribute__((unused)) ElementType el_type,
                                         __attribute__((unused)) GhostType ghost_type = _not_ghost) {}

  /// set the material to steady state (to be implemented for materials that need it)
  virtual void setToSteadyState(__attribute__((unused)) ElementType el_type,
                                __attribute__((unused)) GhostType ghost_type = _not_ghost) {  }

  /// function called to update the internal parameters when the modifiable
  /// parameters are modified
  virtual void updateInternalParameters() {}

public:
  /// extrapolate internal values
  virtual void extrapolateInternal(const ID & id, const Element & element, const Matrix<Real> & points, Matrix<Real> & extrapolated);

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed(const Element & element) const { AKANTU_DEBUG_TO_IMPLEMENT(); }

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed(const Element & element) const { AKANTU_DEBUG_TO_IMPLEMENT(); }

  /// get a material celerity to compute the stable time step (default: is the push wave speed)
  virtual Real getCelerity(const Element & element) const { return getPushWaveSpeed(element); }

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

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// compute the residual for this material
  virtual void updateResidual(GhostType ghost_type = _not_ghost);

  /// assemble the residual for this material
  virtual void assembleResidual(GhostType ghost_type);

  /// Operations before and after solveStep in implicit
  virtual void beforeSolveStep() {}
  virtual void afterSolveStep() {}

  /// save the stress in the previous_stress if needed
  virtual void savePreviousState();

  /// compute the stresses for this material
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);
  virtual void computeAllNonLocalStresses(__attribute__((unused)) GhostType ghost_type = _not_ghost) {};
  virtual void computeAllStressesFromTangentModuli(GhostType ghost_type = _not_ghost);
  virtual void computeAllCauchyStresses(GhostType ghost_type = _not_ghost);

  /// set material to steady state
  void setToSteadyState(GhostType ghost_type = _not_ghost);

  /// compute the stiffness matrix
  virtual void assembleStiffnessMatrix(GhostType ghost_type);

  /// add an element to the local mesh filter
  inline UInt addElement(const ElementType & type,
                         UInt element,
                         const GhostType & ghost_type);

  /// add many elements at once
  void addElements(const Array<Element> & elements_to_add);

  /// remove many element at once
  void removeElements(const Array<Element> & elements_to_remove);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /**
   * interpolate stress on given positions for each element by means
   * of a geometrical interpolation on quadrature points
   */
  void interpolateStress(ElementTypeMapArray<Real> & result,
                         const GhostType ghost_type = _not_ghost);

  /**
   * interpolate stress on given positions for each element by means
   * of a geometrical interpolation on quadrature points and store the
   * results per facet
   */
  void interpolateStressOnFacets(ElementTypeMapArray<Real> & result,
                                 const GhostType ghost_type = _not_ghost);

  /**
   * function to initialize the elemental field interpolation
   * function by inverting the quadrature points' coordinates
   */
  void initElementalFieldInterpolation(const ElementTypeMapArray<Real> & interpolation_points_coordinates);

  /* ------------------------------------------------------------------------ */
  /* Common part                                                              */
  /* ------------------------------------------------------------------------ */
protected:
  /* ------------------------------------------------------------------------ */
  inline UInt getTangentStiffnessVoigtSize(UInt spatial_dimension) const;


  /// compute the potential energy by element
  void computePotentialEnergyByElements();

  /// resize the intenals arrays
  virtual void resizeInternals();

  /* ------------------------------------------------------------------------ */
  /* Finite deformation functions                                             */
  /* This functions area implementing what is described in the paper of Bathe */
  /* et al, in IJNME, Finite Element Formulations For Large Deformation       */
  /* Dynamic Analysis, Vol 9, 353-386, 1975                                   */
  /* ------------------------------------------------------------------------ */
protected:
  /// assemble the residual
  template<UInt dim>
  void assembleResidual(GhostType ghost_type);

  /// Computation of Cauchy stress tensor in the case of finite deformation from
  /// the 2nd Piola-Kirchhoff for a given element type
  template<UInt dim>
  void computeCauchyStress(__attribute__((unused)) ElementType el_type,
                           __attribute__((unused)) GhostType ghost_type = _not_ghost);

  /// Computation the Cauchy stress the 2nd Piola-Kirchhoff and the deformation
  /// gradient
  template<UInt dim >
  inline void computeCauchyStressOnQuad(const Matrix<Real> & F, const Matrix<Real> & S,
                                        Matrix<Real> & cauchy,
                                        const Real & C33 = 1.0) const;

  template<UInt dim>
  void computeAllStressesFromTangentModuli(const ElementType & type,
                                           GhostType ghost_type);

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

  /// Size of the Stress matrix for the case of finite deformation see: Bathe et
  /// al, IJNME, Vol 9, 353-386, 1975
  inline UInt getCauchyStressMatrixSize(UInt spatial_dimension) const;

  /// Sets the stress matrix according to Bathe et al, IJNME, Vol 9, 353-386, 1975
  template<UInt dim>
  inline void setCauchyStressMatrix(const Matrix<Real> & S_t,
                                    Matrix<Real> & sigma);

  /// write the stress tensor in the Voigt notation.
  template<UInt dim>
  inline void setCauchyStressArray(const Matrix<Real> & S_t,
                                   Matrix<Real> & sigma_voight);

  /* ------------------------------------------------------------------------ */
  /* Conversion functions                                                     */
  /* ------------------------------------------------------------------------ */
public:
  template<UInt dim>
  static inline void gradUToF   (const Matrix<Real> & grad_u, Matrix<Real> & F);
  static inline void rightCauchy(const Matrix<Real> & F,      Matrix<Real> & C);
  static inline void leftCauchy (const Matrix<Real> & F,      Matrix<Real> & B);

  template<UInt dim>
  static inline void gradUToEpsilon(const Matrix<Real> & grad_u,
                                    Matrix<Real> & epsilon);
  template<UInt dim>
  static inline void gradUToGreenStrain(const Matrix<Real> & grad_u,
                                        Matrix<Real> & epsilon);

  static inline Real stressToVonMises(const Matrix<Real> & stress);

protected:
  /// converts global element to local element
  inline Element convertToLocalElement(const Element & global_element) const;
  /// converts local element to global element
  inline Element convertToGlobalElement(const Element & local_element) const;

  /// converts global quadrature point to local quadrature point
  inline QuadraturePoint convertToLocalPoint(const QuadraturePoint & global_point) const;
  /// converts local quadrature point to global quadrature point
  inline QuadraturePoint convertToGlobalPoint(const QuadraturePoint & local_point) const;

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
  inline void packElementDataHelper(const ElementTypeMapArray<T> & data_to_pack,
                                    CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const ID & fem_id = ID()) const;

  template<typename T>
  inline void unpackElementDataHelper(ElementTypeMapArray<T> & data_to_unpack,
                                      CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      const ID & fem_id = ID());

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  virtual void onElementsAdded(const Array<Element> & element_list,
                               const NewElementsEvent & event);

  virtual void onElementsRemoved(const Array<Element> & element_list,
                                 const ElementTypeMapArray<UInt> & new_numbering,
                                 const RemovedElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* SolidMechanicsModelEventHandler inherited members                        */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onBeginningSolveStep(const AnalysisMethod & method);
  virtual void onEndSolveStep(const AnalysisMethod & method);
  virtual void onDamageIteration();
  virtual void onDamageUpdate(); 
  virtual void onDump();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(Name,    name,  const std::string &);

  AKANTU_GET_MACRO(Model, *model, const SolidMechanicsModel &)

  AKANTU_GET_MACRO(ID, Memory::getID(), const ID &);
  AKANTU_GET_MACRO(Rho, rho, Real);
  AKANTU_SET_MACRO(Rho, rho, Real);

  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// return the potential energy for the subset of elements contained by the material
  Real getPotentialEnergy();
  /// return the potential energy for the provided element
  Real getPotentialEnergy(ElementType & type, UInt index);

  /// return the energy (identified by id) for the subset of elements contained by the material
  virtual Real getEnergy(std::string energy_id);
  /// return the energy (identified by id) for the provided element
  virtual Real getEnergy(std::string energy_id, ElementType type, UInt index);


  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GradU, gradu, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PotentialEnergy, potential_energy, Real);
  AKANTU_GET_MACRO(GradU, gradu, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(Stress, stress, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(ElementFilter, element_filter, const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(FEEngine, *fem, FEEngine &);

  bool isNonLocal() const { return is_non_local; }

  template <typename T>
  const Array<T> & getArray(const ID & id, const ElementType & type, const GhostType & ghost_type = _not_ghost) const;
  template <typename T>
  Array<T> & getArray(const ID & id, const ElementType & type, const GhostType & ghost_type = _not_ghost);

  template <typename T>
  const InternalField<T> & getInternal(const ID & id) const;
  template <typename T>
  InternalField<T> & getInternal(const ID & id);

  inline bool isInternal(const ID & id, const ElementKind & element_kind) const;
  virtual ElementTypeMap<UInt> getInternalDataPerElem(const ID & id,
                                                      const ElementKind & element_kind,
                                                      const ID & fe_engine_id = "") const;

  bool isFiniteDeformation() const { return finite_deformation; }
  bool isInelasticDeformation() const { return inelastic_deformation; }

  template <typename T>
  inline void setParam(const ID & param, T value);

  template <typename T>
  inline const T & getParam(const ID & param) const;

  virtual void flattenInternal(const std::string & field_id,
                               ElementTypeMapArray<Real> & internal_flat,
                               const GhostType ghost_type = _not_ghost,
                               ElementKind element_kind = _ek_not_defined) const;
protected:
  /// internal variation of the flatten function that is more flexible and can
  /// be used by inherited materials to change some behavior
  virtual void flattenInternalIntern(const std::string & field_id,
                                     ElementTypeMapArray<Real> & internal_flat,
                                     UInt spatial_dimension,
                                     const GhostType ghost_type,
                                     ElementKind element_kind,
                                     const ElementTypeMapArray<UInt> * element_filter = NULL,
                                     const Mesh * mesh = NULL) const;

protected:

  bool isInit() const { return is_init; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// boolean to know if the material has been initialized
  bool is_init;

  std::map<ID, InternalField<Real> *> internal_vectors_real;
  std::map<ID, InternalField<UInt> *> internal_vectors_uint;
  std::map<ID, InternalField<bool> *> internal_vectors_bool;

protected:
  /// Link to the fem object in the model
  FEEngine * fem;

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
  ElementTypeMapArray<UInt> element_filter;

  /// stresses arrays ordered by element types
  InternalField<Real> stress;

  /// eigengrad_u arrays ordered by element types
  InternalField<Real> eigengradu;

  /// grad_u arrays ordered by element types
  InternalField<Real> gradu;

  /// Green Lagrange strain (Finite deformation)
  InternalField<Real> green_strain;

  /// Second Piola-Kirchhoff stress tensor arrays ordered by element types (Finite deformation)
  InternalField<Real> piola_kirchhoff_2;

  /// potential energy by element
  InternalField<Real> potential_energy;

  /// tell if using in non local mode or not
  bool is_non_local;

  /// tell if the material need the previous stress state
  bool use_previous_stress;

  /// tell if the material need the previous strain state
  bool use_previous_gradu;

  /// elemental field interpolation coordinates
  InternalField<Real> interpolation_inverse_coordinates;

  /// elemental field interpolation points
  InternalField<Real> interpolation_points_matrices;

  /// vector that contains the names of all the internals that need to
  /// be transferred when material interfaces move
  std::vector<ID> internals_to_transfer;
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
/// This can be used to automatically write the loop on quadrature points in
/// functions such as computeStress. This macro in addition to write the loop
/// provides two tensors (matrices) sigma and grad_u
#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type) \
  Array<Real>::matrix_iterator gradu_it =                               \
    this->gradu(el_type, ghost_type).begin(this->spatial_dimension,     \
                                          this->spatial_dimension);     \
  Array<Real>::matrix_iterator gradu_end =                              \
    this->gradu(el_type, ghost_type).end(this->spatial_dimension,       \
                                        this->spatial_dimension);       \
                                                                        \
  this->stress(el_type,                                                 \
               ghost_type).resize(this->gradu(el_type,                  \
                                              ghost_type).getSize());   \
                                                                        \
  Array<Real>::iterator< Matrix<Real> > stress_it =                     \
    this->stress(el_type, ghost_type).begin(this->spatial_dimension,    \
                                            this->spatial_dimension);   \
                                                                        \
  if(this->isFiniteDeformation()){                                      \
    this->piola_kirchhoff_2(el_type,                                    \
                            ghost_type).resize(this->gradu(el_type,     \
                                                           ghost_type).getSize()); \
    stress_it =                                                         \
      this->piola_kirchhoff_2(el_type,                                  \
                              ghost_type).begin(this->spatial_dimension, \
                                                this->spatial_dimension); \
  }                                                                     \
                                                                        \
  for(;gradu_it != gradu_end; ++gradu_it, ++stress_it) {                \
    Matrix<Real> & __attribute__((unused)) grad_u = *gradu_it;          \
    Matrix<Real> & __attribute__((unused)) sigma  = *stress_it

#define MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END                       \
  }                                                                     \


/// This can be used to automatically write the loop on quadrature points in
/// functions such as computeTangentModuli. This macro in addition to write the
/// loop provides two tensors (matrices) sigma_tensor, grad_u, and a matrix
/// where the elemental tangent moduli should be stored in Voigt Notation
#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_mat)       \
  Array<Real>::matrix_iterator gradu_it =                               \
    this->gradu(el_type, ghost_type).begin(this->spatial_dimension,     \
                                           this->spatial_dimension);    \
  Array<Real>::matrix_iterator gradu_end =                              \
    this->gradu(el_type, ghost_type).end(this->spatial_dimension,       \
                                         this->spatial_dimension);      \
  Array<Real>::matrix_iterator sigma_it =                               \
    this->stress(el_type, ghost_type).begin(this->spatial_dimension,    \
                                            this->spatial_dimension);   \
                                                                        \
  tangent_mat.resize(this->gradu(el_type, ghost_type).getSize());       \
                                                                        \
  UInt tangent_size =                                                   \
    this->getTangentStiffnessVoigtSize(this->spatial_dimension);        \
  Array<Real>::matrix_iterator tangent_it =                             \
    tangent_mat.begin(tangent_size,                                     \
                      tangent_size);                                    \
                                                                        \
  for(;gradu_it != gradu_end; ++gradu_it, ++sigma_it, ++tangent_it) {   \
    Matrix<Real> & __attribute__((unused)) grad_u  = *gradu_it;         \
    Matrix<Real> & __attribute__((unused)) sigma_tensor = *sigma_it;    \
    Matrix<Real> & tangent = *tangent_it


#define MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END                      \
  }                                                                     \

/* -------------------------------------------------------------------------- */
#define INSTANTIATE_MATERIAL(mat_name)                  \
  template class mat_name<1>;                           \
  template class mat_name<2>;                           \
  template class mat_name<3>

#endif /* __AKANTU_MATERIAL_HH__ */
