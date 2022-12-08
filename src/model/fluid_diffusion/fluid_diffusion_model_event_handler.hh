/**
 * @file   fluid_diffusion_model_event_handler.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Mon Aug 19 2019
 * @date last modification:
 *
 * @brief  EventHandler implementation for FluidDiffusionEvents
 * modification of the SolidMecanicsModelEventHandler
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

#ifndef __AKANTU_FLUID_DIFFUSION_MODEL_EVENT_HANDLER_HH__
#define __AKANTU_FLUID_DIFFUSION_MODEL_EVENT_HANDLER_HH__

namespace akantu {

/// akantu::FluidDiffusionModelEvent is the base event for model
namespace FluidDiffusionModelEvent {}

/// akantu::SolidMechanicsModelEvent
class FluidDiffusionModelEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual ~FluidDiffusionModelEventHandler() = default;

  template <class EventHandler> friend class EventHandlerManager;

};

} // namespace akantu

#endif /* __AKANTU_FLUID_DIFFUSION_MODEL_EVENT_HANDLER_HH__ */
