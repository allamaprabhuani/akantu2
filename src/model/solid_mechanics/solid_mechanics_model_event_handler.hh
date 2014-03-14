/**
 * @file   solid_mechanics_model_event_handler.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Wed Mar 12 16:48:21 2014
 *
 * @brief
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

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_EVENT_HANDLER_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_EVENT_HANDLER_HH__

__BEGIN_AKANTU__

namespace SolidMechanicsModelEvent {
  struct BeforeSolveStepEvent {
    BeforeSolveStepEvent(AnalysisMethod & method) : method(method) {}
    AnalysisMethod method;
  };
  struct AfterSolveStepEvent {
    AfterSolveStepEvent(AnalysisMethod & method) : method(method) {}
    AnalysisMethod method;
  };

}


class SolidMechanicsModelEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual ~SolidMechanicsModelEventHandler() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:

  inline void sendEvent(const SolidMechanicsModelEvent::BeforeSolveStepEvent & event) {
    onBeginningSolveStep(event.method);
  }
  inline void sendEvent(const SolidMechanicsModelEvent::AfterSolveStepEvent & event) {
    onEndSolveStep(event.method);
  }

  template<class EventHandler>
  friend class EventHandlerManager;

  /* ------------------------------------------------------------------------ */
  /* Interface                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onBeginningSolveStep(__attribute__((unused)) const AnalysisMethod & method) {}
  virtual void onEndSolveStep(__attribute__((unused)) const AnalysisMethod & method) {}
};


__END_AKANTU__

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_EVENT_HANDLER_HH__ */
