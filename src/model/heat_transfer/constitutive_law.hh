/**
 * @file   constitutive_law.hh
 *
 * @author Mohit Pundir <mohit.pundir@ethz.ch>
 *
 * @date creation: Mon May 2 2022
 * @date last modification: Wed may 2 2022
 *
 * @brief  Mother class for all constitutive laws
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_factory.hh"
#include "data_accessor.hh"
#include "parsable.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
#include "internal_field.hh"
#include "random_internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONSTITUTIVE_LAW_HH__
#define __AKANTU_CONSTITUTIVE_LAW_HH__

namespace akantu{
  class Model;
  class PoissonModel;
  class ConstitutiveLaw;
} // namespace akantu


namespace akantu{

template <typename T>
using InternalConstitutiveLaw = InternalFieldTmpl<ConstitutiveLaw, T>;

using ConstitutiveLawFactory =
    Factory<ConstitutiveLaw, ID, const ID &, PoissonModel &, const ID &>;

class ConstitutiveLaw : public DataAccessor<Element>, public Parsable {

  
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ConstitutiveLaw(const ConstitutiveLaw & law) = delete;
  ConstitutiveLaw & operator=(const ConstitutiveLaw & law) = delete;

  /// Initialize constitutive law with defaults
  ConstitutiveLaw(PoissonModel & model, const ID & id = "");
  
  /// Initialize constitutive law with custom mesh & fe_engine
  ConstitutiveLaw(PoissonModel & model, UInt dim, const Mesh & mesh,
		  FEEngine & fe_engine, const ID & id = "");

  /// Destructor
  ~ConstitutiveLaw() override;


protected:
  void initialize();

  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename T>
  void registerInternal(InternalConstitutiveLaw<T> & /*vect*/) {
    AKANTU_TO_IMPLEMENT();
  }

  template <typename T>
  void unregisterInternal(InternalConstitutiveLaw<T> & /*vect*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// initialize the constitutive law computed parameter
  virtual void initConstitutiveLaw();

  ///
  virtual void beforeSolveStep();

  ///
  virtual void afterSolveStep(bool converged = true);

  /// compute the fluxes for this constitutive law
  virtual void computeAllFluxes(GhostType ghost_type = _not_ghost);
  
    /// assemble the internal dof rate for this constitutive law
  virtual void assembleInternalDofRate(GhostType ghost_type);

  /// compute the stiffness matrix
  virtual void assembleStiffnessMatrix(GhostType ghost_type);
   
  /// add an element to the local mesh filter
  inline UInt addElement(const ElementType & type, UInt element,
                         const GhostType & ghost_type);
  inline UInt addElement(const Element & element);

  /// get a constitutive law celerity to compute the stable time step
  virtual Real getCelerity() const {
    AKANTU_TO_IMPLEMENT();
  }

  /// get a effective capacity to assmble capacity matrix
  virtual Real getEffectiveCapacity() const {
    AKANTU_TO_IMPLEMENT();
  }

  
  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;
 

protected:
  /// resize the internals arrrays
  virtual void resizeInternals();

  /// function called to updatet the internal parameters when the
  /// modifiable parameters are modified
  virtual void updateInternalParameters() {
    AKANTU_TO_IMPLEMENT();
  }

  
  /* ------------------------------------------------------------------------ */
  /* Function that constitutive laws can/should reimplement                           */
  /* ------------------------------------------------------------------------ */
protected:
  
  /// compute flux
  virtual void computeFlux(ElementType /* el_type */,
                             GhostType /* ghost_type */ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }


  /// compute tangent modulii
  virtual void computeTangentModuli(ElementType /*el_type*/,
                                    Array<Real> & /*tangent_matrix*/,
                                    GhostType /*ghost_type*/ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }
  
  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  template <typename T>
  inline void packElementDataHelper(const ElementTypeMapArray<T> & data_to_pack,
                                    CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const ID & fem_id = ID()) const;

  template <typename T>
  inline void unpackElementDataHelper(ElementTypeMapArray<T> & data_to_unpack,
                                      CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      const ID & fem_id = ID());

  

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Name, name, const std::string &);
  AKANTU_SET_MACRO(Name, name, const std::string &);

  AKANTU_GET_MACRO(Model, model, const PoissonModel &)

  AKANTU_GET_MACRO(ID, id, const ID &);
  
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(FluxDof, flux_dof, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(GradientDof, gradient_dof, Real);
  
  AKANTU_GET_MACRO(GradientDof, gradient_dof, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(FluxDof, flux_dof, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(ElementFilter, element_filter,
                   const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(FEEngine, fem, FEEngine &);

  
  
  template <typename T>
  const InternalConstitutiveLaw<T> & getInternal(const ID & id) const;
  template <typename T> InternalConstitutiveLaw<T> & getInternal(const ID & id);

  template <typename T>
  inline bool isInternal(const ID & id, const ElementKind & element_kind) const;

  template <typename T> inline void setParam(const ID & param, T value);
  inline const Parameter & getParam(const ID & param) const;

  template <typename T>
  void flattenInternal(const std::string & field_id,
                       ElementTypeMapArray<T> & internal_flat,
                       GhostType ghost_type = _not_ghost,
                       ElementKind element_kind = _ek_not_defined) const;

  
  bool hasMatrixChanged(const ID & id) {
    if (id == "K") {
      return hasStiffnessMatrixChanged();
    }

    return true;
  }


  MatrixType getMatrixType(const ID & id) {
    if (id == "K" or id == "M") {
      return _symmetric;
    }

    return _mt_not_defined;
  }

  
  /// static method to reteive the material factory
  static ConstitutiveLawFactory & getFactory();


  /// specify if the matrix need to be recomputed for this
  /// constitutive law
  virtual bool hasStiffnessMatrixChanged() { return true; }


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// boolean to know if the material has been initialized
  bool is_init;

  std::map<ID, InternalConstitutiveLaw<Real> *> internal_vectors_real;
  std::map<ID, InternalConstitutiveLaw<UInt> *> internal_vectors_uint;
  std::map<ID, InternalConstitutiveLaw<bool> *> internal_vectors_bool;

protected:
  ID id;

  /// Link to the fem object in the model
  FEEngine & fem;

  /// phasefield name
  std::string name;

  /// The model to whch the constitutive Law belong
  PoissonModel & model;

  /// spatial dimension
  UInt spatial_dimension;

  /// list of element handled by the constitutive law
  ElementTypeMapArray<UInt> element_filter;

  /// fluxes arrays ordered by element types
  InternalConstitutiveLaw<Real> flux_dof;

  /// gradient of dof arrays ordered by element types
  InternalConstitutiveLaw<Real> gradient_dof;
  
  /// vector that contains the names of all the internals that need to
  /// be transferred when constitutive law interfaces move
  std::vector<ID> internals_to_transfer;
  
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const ConstitutiveLaw & _this) {
  _this.printself(stream);
  return stream;
}
  

} // namespace akantu

  
#include "constitutive_law_inline_impl.hh"


#include "internal_field_tmpl.hh"
#include "random_internal_field_tmpl.hh"

/* -------------------------------------------------------------------------- */
#define CONSTITUTIVE_LAW_DEFAULT_ALLOCATOR(id, law_name)		\
  [](const ID &, PoissonModel & model,					\
     const ID & id) -> std::unique_ptr<ConstitutiveLaw> {		\
    return std::make_unique<law_name>(model, id);			\
  }

#define INSTANTIATE_CONSTITUTIVE_LAW(id, law_name)			\
  static bool constitutive_law_is_alocated_##id [[gnu::unused]] =	\
    ConstitutiveLawFactory::getInstance().registerAllocator(		\
	 #id, CONSTITUTIVE_LAW_DEFAULT_ALLOCATOR(id, law_name))

#endif
