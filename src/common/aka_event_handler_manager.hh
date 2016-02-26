/**
 * @file   aka_event_handler_manager.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Dec 16 2015
 *
 * @brief  Base of Event Handler classes
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

#ifndef __AKANTU_AKA_EVENT_HANDLER_MANAGER_HH__
#define __AKANTU_AKA_EVENT_HANDLER_MANAGER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <list>
#include <algorithm>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template <class EventHandler> class EventHandlerManager {
private:
  typedef std::pair<UInt, EventHandler *> priority_value;
  typedef std::list<priority_value> priority_list;
  struct KeyComp {
    bool operator()(const priority_value & a, const priority_value & b) const {
      return (a.first < b.first);
    }
    bool operator()(const priority_value & a, UInt b) const {
      return (a.first < b);
    }

  };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual ~EventHandlerManager(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register a new EventHandler to the Manager. The register object
  /// will then be informed about the events the manager observes.
  void registerEventHandler(EventHandler & event_handler, UInt priority = 100) {
    typename priority_list::iterator it =
        this->searchEventHandler(event_handler);

    if (it != this->event_handlers.end()) {
      AKANTU_EXCEPTION("This event handler was already registered");
    }

    typename priority_list::iterator pos =
        std::lower_bound(this->event_handlers.begin(),
                         this->event_handlers.end(), priority, KeyComp());

    this->event_handlers.insert(pos, std::make_pair(priority, &event_handler));
  }

  /// unregister a EventHandler object. This object will not be
  /// notified anymore about the events this manager observes.
  void unregisterEventHandler(EventHandler & event_handler) {
    typename priority_list::iterator it =
        this->searchEventHandler(event_handler);

    if (it == this->event_handlers.end()) {
      AKANTU_EXCEPTION("This event handler is not registered");
    }

    this->event_handlers.erase(it);
  }

  /// Notify all the registered EventHandlers about the event that just occured.
  template <class Event> void sendEvent(const Event & event) {
    typename priority_list::iterator it = event_handlers.begin();
    typename priority_list::iterator end = event_handlers.end();
    for (; it != end; ++it)
      it->second->sendEvent(event);
  }

private:
  typename priority_list::iterator searchEventHandler(EventHandler & handler) {
    typename priority_list::iterator it = this->event_handlers.begin();
    typename priority_list::iterator end = this->event_handlers.end();

    for (; it != end && it->second != &handler; ++it)
      ;

    return it;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// list of the event handlers
  priority_list event_handlers;
};

__END_AKANTU__

#endif /* __AKANTU_AKA_EVENT_HANDLER_MANAGER_HH__ */
