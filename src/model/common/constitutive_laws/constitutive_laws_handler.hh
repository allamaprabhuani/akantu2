/**
 * @file   constitutive_laws_handler.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation  ven mar 26 2021
 *
 * @brief A handler for multiple consistutive laws
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
#include "constitutive_law.hh"
#include "constitutive_law_selector.hh"
#include "mesh.hh"
#include "mesh_events.hh"
#include "non_local_manager_callback.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONSTITUTIVE_LAWS_HANDLER_HH
#define AKANTU_CONSTITUTIVE_LAWS_HANDLER_HH

namespace akantu {
class NonLocalManager;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

template <class ConstitutiveLawType>
class ConstitutiveLawsHandler : public DataAccessor<Element>,
                                public NonLocalManagerCallback {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ConstitutiveLawsHandler(const Mesh & mesh, UInt spatial_dimension,
                          const ID & parent_id)
      : constitutive_law_index("constitutive_law index", parent_id),
        constitutive_law_local_numbering("constitutive_law local numbering",
                                         parent_id) {
    constitutive_law_selector =
        std::make_shared<DefaultConstitutiveLawSelector>(
            constitutive_law_index);

    constitutive_law_index.initialize(mesh, _element_kind = _ek_not_defined,
                                      _default_value = UInt(-1),
                                      _with_nb_element = true);
    constitutive_law_local_numbering.initialize(
        mesh, _element_kind = _ek_not_defined, _with_nb_element = true);
  }

  virtual ~ConstitutiveLawsHandler();

  class NewConstitutiveLawElementsEvent : public NewElementsEvent {
  public:
    AKANTU_GET_MACRO_NOT_CONST(ConstitutiveLawList, constitutive_law,
                               Array<UInt> &);
    AKANTU_GET_MACRO(ConstitutiveLawList, constitutive_law,
                     const Array<UInt> &);

  protected:
    Array<UInt> constitutive_law;
  };

  /* ------------------------------------------------------------------------ */
  /* ConstitutiveLaws                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// register an empty constitutive_law of a given type
  auto & registerNewConstitutiveLaw(const ID & cl_name, const ID & cl_type,
                                    const ID & opt_param);

  /// reassigns constitutive_laws depending on the constitutive_law selector
  virtual void reassignConstitutiveLaw();

  template <class Func> void for_each_constitutive_law(Func && func) {
    for (auto && constitutive_law : constitutive_laws) {
      std::forward<Func>(func)(
          aka::as_type<ConstitutiveLawType>(*constitutive_laws));
    }
  }

protected:
  /// register a constitutive_law in the dynamic database
  auto & registerNewConstitutiveLaw(const ParserSection & cl_section);

  /// read the constitutive_law files to instantiate all the constitutive_laws
  void instantiateConstitutiveLaws();

  /// set the element_id_by_constitutive_law and add the elements to the good
  /// constitutive_laws
  virtual void assignConstitutiveLawToElements(
      const ElementTypeMapArray<UInt> * filter = nullptr);

protected:
  void splitElementByConstitutiveLaw(
      const Array<Element> & elements,
      std::vector<Array<Element>> & elements_per_cl) const;

  template <typename Operation>
  void splitByConstitutiveLaw(const Array<Element> & elements,
                              Operation && op) const;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface (kept for convenience) and dumper relative functions  */
  /* ------------------------------------------------------------------------ */
public:
  //! decide wether a field is a constitutive_law internal or not
  bool isInternal(const std::string & field_name, ElementKind element_kind);

  //! give the amount of data per element
  virtual ElementTypeMap<UInt>
  getInternalDataPerElem(const std::string & field_name, ElementKind kind);

  //! flatten a given constitutive_law internal field
  ElementTypeMapArray<Real> &
  flattenInternal(const std::string & field_name, ElementKind kind,
                  GhostType ghost_type = _not_ghost);

  //! flatten all the registered constitutive_law internals
  void flattenAllRegisteredInternals(ElementKind kind);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get an iterable on the constitutive_laws
  inline decltype(auto) getConstitutiveLaws();

  /// get an iterable on the constitutive_laws
  inline decltype(auto) getConstitutiveLaws() const;

  /// get a particular constitutive_law (by numerical constitutive_law index)
  inline auto & getConstitutiveLaw(UInt cl_index);

  /// get a particular constitutive_law (by numerical constitutive_law index)
  inline const auto & getConstitutiveLaw(UInt cl_index) const;

  /// get a particular constitutive_law (by constitutive_law name)
  inline auto & getConstitutiveLaw(const std::string & name);

  /// get a particular constitutive_law (by constitutive_law name)
  inline const auto & getConstitutiveLaw(const std::string & name) const;

  /// get a particular constitutive_law id from is name
  inline UInt getConstitutiveLawIndex(const std::string & name) const;

  /// give the number of constitutive_laws
  inline UInt getNbConstitutiveLaws() const { return constitutive_laws.size(); }

  /// give the constitutive_law internal index from its id
  Int getInternalIndexFromID(const ID & id) const;

  AKANTU_GET_MACRO(ConstitutiveLawByElement, constitutive_law_index,
                   const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(ConstitutiveLawLocalNumbering,
                   constitutive_law_local_numbering,
                   const ElementTypeMapArray<UInt> &);

  /// vectors containing local constitutive_law element index for each global
  /// element index
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConstitutiveLawByElement,
                                         constitutive_law_index, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConstitutiveLawLocalNumbering,
                                         constitutive_law_local_numbering,
                                         UInt);

  AKANTU_GET_MACRO_NOT_CONST(ConstitutiveLawSelector,
                             *constitutive_law_selector,
                             ConstitutiveLawSelector &);

  /// Access the non_local_manager interface
  AKANTU_GET_MACRO(NonLocalManager, *non_local_manager, NonLocalManager &);

  void setConstitutiveLawSelector(
      std::shared_ptr<ConstitutiveLawSelector> constitutive_law_selector) {
    this->constitutive_law_selector = std::move(constitutive_law_selector);
  }

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  UInt getNbData(const Array<Element> & elements,
                 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  void onElementsAdded(const Array<Element> & element_list,
                       const NewElementsEvent & event) override;
  void onElementsRemoved(const Array<Element> & element_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const RemovedElementsEvent & event) override;

  /* ------------------------------------------------------------------------ */
  /* NonLocalManager inherited members                                        */
  /* ------------------------------------------------------------------------ */
protected:
  void initializeNonLocal() override;

  void updateDataForNonLocalCriterion(ElementTypeMapReal & criterion) override;

  void computeNonLocalStresses(GhostType ghost_type) override;

  void insertIntegrationPointsInNeighborhoods(GhostType ghost_type) override;

  /// update the values of the non local internal
  void updateLocalInternal(ElementTypeMapReal & internal_flat,
                           GhostType ghost_type, ElementKind kind) override;

  /// copy the results of the averaging in the constitutive_laws
  void updateNonLocalInternal(ElementTypeMapReal & internal_flat,
                              GhostType ghost_type, ElementKind kind) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// id of the derived class
  const ID & parent_id;

  /// underlying mesh
  const Mesh & _mesh;

  /// spatial_dimension of the derived model
  UInt _spatial_dimension{_all_dimensions};

  /// mapping between constitutive_law name and constitutive_law internal id
  std::map<std::string, UInt> constitutive_laws_names_to_id;

  /// Arrays containing the constitutive_law index for each element
  ElementTypeMapArray<UInt> constitutive_law_index;

  /// Arrays containing the position in the element filter of the
  /// constitutive_law (constitutive_law's local numbering)
  ElementTypeMapArray<UInt> constitutive_law_local_numbering;

  /// list of used constitutive_laws
  std::vector<std::unique_ptr<ConstitutiveLawType>> constitutive_laws;

  /// class defining of to choose a constitutive_law
  std::shared_ptr<ConstitutiveLawSelector> constitutive_law_selector;

  using flatten_internal_map =
      std::map<std::pair<std::string, ElementKind>,
               std::unique_ptr<ElementTypeMapArray<Real>>>;

  /// map a registered internals to be flattened for dump purposes
  flatten_internal_map registered_internals;

  /// tells if the constitutive_law are instantiated
  bool are_constitutive_laws_instantiated{false};

  /// non local manager
  std::unique_ptr<NonLocalManager> non_local_manager;

  template <class ConstitutiveLawsHandler> friend class ConstitutiveLaw;
};
  
} // namespace akantu

#include "constitutive_laws_handler_tmpl.hh"
#endif /* AKANTU_CONSTITUTIVE_LAWS_HANDLER_HH */
