/**
 * @file   contact_mechanics_model.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Wed Feb 21 2018
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
#include "model.hh"
#include "contact_detector.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_MECHANICS_MODEL_HH__
#define __AKANTU_CONTACT_MECHANICS_MODEL_HH__

namespace akantu {
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class ContactMechanicsModel
  : public Model {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  
  
public:
  ContactMechanicsModel(
      Mesh & mesh, UInt spatial_dimension = _all_dimensions,
      const ID & id = "contact_mechanics_model", const MemoryID & memory_id = 0,
      const ModelType model_type = ModelType::_contact_mechanics_model);

  ~ContactMechanicsModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:

  /// initialize the modelType
  void initModel() override;


  friend class ContactDetector;
  
private:

  void computeNormals();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// contact detection
  std::unique_ptr<ContactDetector> detector;
  
};




} // namespace akantu



#endif
