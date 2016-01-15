/**
 * @file   solid_mechanics_model_event_handler.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Dec 18 2015
 *
 * @brief  EventHandler implementation for SolidMechanicsEvents
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

///akantu::SolidMechanicsModelEvent is the base event for model
namespace SolidMechanicsModelEvent {
  struct BeforeSolveStepEvent {
    BeforeSolveStepEvent(AnalysisMethod & method) : method(method) {}
    AnalysisMethod method;
  };
  struct AfterSolveStepEvent {
    AfterSolveStepEvent(AnalysisMethod & method) : method(method) {}
    AnalysisMethod method;
  };
  struct BeforeDumpEvent { BeforeDumpEvent() {} };
  struct BeginningOfDamageIterationEvent { BeginningOfDamageIterationEvent() {} };
  struct AfterDamageEvent { AfterDamageEvent() {} };
}

/// akantu::SolidMechanicsModelEvent
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
  ///Send what is before the solve step to the beginning of solve step through EventManager
  inline void sendEvent(const SolidMechanicsModelEvent::BeforeSolveStepEvent & event) {
    onBeginningSolveStep(event.method);
  }
  ///Send what is after the solve step to the end of solve step through EventManager
  inline void sendEvent(const SolidMechanicsModelEvent::AfterSolveStepEvent & event) {
    onEndSolveStep(event.method);
  }
  ///Send what is before dump to current dump through EventManager
  inline void sendEvent(const SolidMechanicsModelEvent::BeforeDumpEvent & event) {
    onDump();
  }
  ///Send what is at the beginning of damage iteration to Damage iteration through EventManager
  inline void sendEvent(const SolidMechanicsModelEvent::BeginningOfDamageIterationEvent & event) {
    onDamageIteration();
  } 
  ///Send what is after damage for the damage update through EventManager
  inline void sendEvent(const SolidMechanicsModelEvent::AfterDamageEvent & event) {
    onDamageUpdate();
  }

  template<class EventHandler>
  friend class EventHandlerManager;

  /* ------------------------------------------------------------------------ */
  /* Interface                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// function to implement to react on akantu::BeforeSolveStepEvent
  virtual void onBeginningSolveStep(__attribute__((unused)) const AnalysisMethod & method) {}
  /// function to implement to react on akantu::AfterSolveStepEvent
  virtual void onEndSolveStep(__attribute__((unused)) const AnalysisMethod & method) {}
  /// function to implement to react on akantu::BeforeDumpEvent
  virtual void onDump() {}
  /// function to implement to react on akantu::BeginningOfDamageIterationEvent
  virtual void onDamageIteration() {}
  /// function to implement to react on akantu::AfterDamageEvent
  virtual void onDamageUpdate() {}
};


__END_AKANTU__

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_EVENT_HANDLER_HH__ */
