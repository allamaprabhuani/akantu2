/**
 * @file   contact_mechanics_model.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Sep 10 2018
 * @date last modification: Mon Sep 10 2018
 *
 * @brief  Model of Contact Mechanics
 *
 * @section LICENSE
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
#include "aka_named_arguments.hh"
#include "model.hh"
#include "contact_resolution.hh"
#include "contact_detection.hh"
#include "contact_element.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_CONTACT_MECHANICS_MODEL_HH__
#define __AKANTU_CONTACT_MECHANICS_MODEL_HH__

namespace akantu {

template<Model model>
class ContactMechanicsModel : public Memory {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ContactMechanicsModel(const ID & id = "contact_mechanics_model",
			const MemoryID & memory_id = 0);

  ~ContactMechanicsModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize completely the model
  void initFullImpl(const ContactModelOptions & options);

  ///
  void initModel();

  /// function to print the content of the class
  void printself(std::ostream &, int = 0) const;
  
public:
#ifndef SWIG
  template <typename... pack>
  std::enable_if_t<are_named_argument<pack...>::value>
  initFull(pack &&... _pack) {
    this->initFullImpl(ContactModelOptions{
	use_named_args, std::forward<decltype(_pack)>(_pack)...});
  }

  template <typename... pack>
  std::enable_if_t<not are_named_argument<pack...>::value>
  initFull(pack &&... _pack) {
    this->initFullImpl(std::forward<decltype(_pack)>(_pack)...);
  }
#endif
  
  /// initialize a new resolution if needed
  void initNewResolution(const ContactResolutionMethod & resolution_method);

  /// create contact elements
  void createContactElements();

  /// update contact elements after resolution
  void updateContactElements();
  
  ///
  bool solve();
  

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
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

protected:
    
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// mesh inherited from class Model
  Mesh & mesh;

  /// 
  std::vector<ContactElement> element;
   
  /// resolution method check the list in akantu::ContactresolutionMethod
  ContactResolutionMethod resolution_method;
  
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const ContactMechanicsModel & _this) {
  _this.printself(stream);
  return stream;
}

  
namespace {
  DECLARE_NAMED_ARGUMENT(resolution_method);
}
  
struct ContactModelOptions {
  explicit ContactModelOptions(ContactResolutionMethod resolution_method =
			       ContactResolutionMethod::_penalty)  
    : resolution_method(resolution_method) {}

  template<typename... pack>    
  ContactModelOptions(use_named_args_t, pack &&... _pack)
    : ContactModelOptions(OPTIONAL_NAMED_ARGS(resolution_method,
					      ContactResolutionMethod::_penalty)) {}

  virtual ~ContactModelOptions() = default;

  ContactResolutionMethod resolution_method;
};
  

}

#endif /* __AKANTU_CONTACT_MECHANICS_MODEL_HH__
