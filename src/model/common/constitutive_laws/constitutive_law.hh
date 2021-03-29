/**
 * @file   constitutive_law.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Mother class for all constitutive laws
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
#include "mesh_events.hh"
#include "parsable.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
#include "internal_field.hh"
#include "random_internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONSTITUTIVE_LAW_HH__
#define __AKANTU_CONSTITUTIVE_LAW_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {

class ConstitutiveLawInternalHandler {
public:
  ConstitutiveLawInternalHandler(const ID & id) : id(id) {}

  template <typename T>
  inline void registerInternal(std::shared_ptr<InternalField<T>> & vect);

  inline void unregisterInternal(const ID & id);

  /// resize the internals arrrays
  virtual void resizeInternals();

public:
  /// save the internals in the previous_state if needed
  virtual void savePreviousState();

  /// restore the internals from previous_state if needed
  virtual void restorePreviousState();

  template <typename T>
  const InternalField<T> & getInternal(const ID & id) const;

  template <typename T> InternalField<T> & getInternal(const ID & id);

  template <typename T>
  inline bool isInternal(const ID & id, const ElementKind & element_kind) const;

  template <typename T>
  const Array<T> & getArray(const ID & id, ElementType type,
                            GhostType ghost_type = _not_ghost) const;
  template <typename T>
  Array<T> & getArray(const ID & id, ElementType type,
                      GhostType ghost_type = _not_ghost);

public:
  virtual FEEngine & getFEEngine() = 0;
  virtual UInt getSpatialDimension() = 0;
  AKANTU_GET_MACRO(Name, name, const std::string &);
  AKANTU_GET_MACRO(ID, id, const ID &);

private:
  std::map<ID, std::shared_ptr<InternalFieldBase>> internal_vectors;

protected:
  ID id;

  /// constitutive law name
  std::string name;
};

template <class ConstitutiveLawsHandler_>
class ConstitutiveLaw : public ConstitutiveLawInternalHandler,
                        public DataAccessor<Element>,
                        public MeshEventHandler,
                        public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using ConstitutiveLawsHandler = ConstitutiveLawsHandler_;

  ConstitutiveLaw(const ConstitutiveLaw & law) = delete;

  ConstitutiveLaw & operator=(const ConstitutiveLaw & law) = delete;

  /// Initialize constitutive law with defaults
  ConstitutiveLaw(ConstitutiveLawsHandler & handler, const ID & id = "",
                  ElementKind element_kind = _ek_regular);

  /// Destructor
  ~ConstitutiveLaw() override;

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the constitutive law computed parameter
  virtual void initConstitutiveLaw();

  /// add an element to the local mesh filter
  inline UInt addElement(ElementType type, UInt element, GhostType ghost_type);
  inline UInt addElement(const Element & element);

  /// add many elements at once
  void addElements(const Array<Element> & elements_to_add);

  /// remove many element at once
  void removeElements(const Array<Element> & elements_to_remove);

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  /// function called to update the internal parameters when the
  /// modifiable parameters are modified
  virtual void updateInternalParameters() {}

  /// converts global element to local element
  inline Element convertToLocalElement(const Element & global_element) const;
  /// converts local element to global element
  inline Element convertToGlobalElement(const Element & local_element) const;

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
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
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  virtual void onNodesAdded(const Array<UInt> & /*unused*/,
                            const NewNodesEvent & /*unused*/) override{};
  virtual void onNodesRemoved(const Array<UInt> & /*unused*/,
                              const Array<UInt> & /*unused*/,
                              const RemovedNodesEvent & /*unused*/) override{};
  virtual void
  onElementsChanged(const Array<Element> & /*unused*/,
                    const Array<Element> & /*unused*/,
                    const ElementTypeMapArray<UInt> & /*unused*/,
                    const ChangedElementsEvent & /*unused*/) override{};

  virtual void onElementsAdded(const Array<Element> & /*unused*/,
                               const NewElementsEvent & /*unused*/);

  virtual void
  onElementsRemoved(const Array<Element> & element_list,
                    const ElementTypeMapArray<UInt> & new_numbering,
                    const RemovedElementsEvent & event);

public:
  template <typename T> inline void setParam(const ID & param, T value);
  inline const Parameter & getParam(const ID & param) const;

  template <typename T>
  void flattenInternal(const std::string & field_id,
                       ElementTypeMapArray<T> & internal_flat,
                       const GhostType ghost_type = _not_ghost,
                       ElementKind element_kind = _ek_not_defined) const;

  /* ------------------------------------------------------------------------
   */
  /* Accessors */
  /* ------------------------------------------------------------------------
   */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ElementFilter, element_filter, UInt);

  template <typename T>
  ElementTypeMap<UInt> getInternalDataPerElem(const ID & id,
                                              ElementKind element_kind) const;

protected:
  bool isInit() const { return is_init; }

  /* ------------------------------------------------------------------------
   */
  /* Class Members */
  /* ------------------------------------------------------------------------
   */
private:
  /// boolean to know if the constitutive law has been initialized
  bool is_init{false};

protected:
  /// list of element handled by the constitutive law
  ElementTypeMapArray<UInt> element_filter;

  // Constitutive law handler for which this is a constitutive law
  ConstitutiveLawsHandler & handler;
};

} // namespace akantu

#include "constitutive_law_tmpl.hh"

#endif
