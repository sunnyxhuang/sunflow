
#include <iomanip>
#include <map>
#include "events.h"
#include "trafficGen.h"
#include "scheduler.h" 
#include "util.h"

using namespace std;

ostream &operator<<(std::ostream &out, const EventType value) {
  static map<EventType, string> strings;
  if (strings.size() == 0) {
#define INSERT_ELEMENT(p) strings[p] = #p
    INSERT_ELEMENT(NO_EVENT);
    INSERT_ELEMENT(SUB_JOB);
    INSERT_ELEMENT(WORKER_FINISH);
    INSERT_ELEMENT(MSG_ADD_FLOWS);
    INSERT_ELEMENT(MSG_ADD_COFLOWS);
    INSERT_ELEMENT(ALARM_TRAFFIC);
    INSERT_ELEMENT(FLOW_FINISH);
    INSERT_ELEMENT(COFLOW_ARRIVE);
    INSERT_ELEMENT(FLOW_ARRIVE);
    INSERT_ELEMENT(RESCHEDULE);
    INSERT_ELEMENT(ORDER_CIRCUIT);
    INSERT_ELEMENT(APPLY_CIRCUIT);
    INSERT_ELEMENT(APPLY_NEW_SCHEDULE);
    INSERT_ELEMENT(SCHEDULE_END);
    INSERT_ELEMENT(MSG_TRAFFIC_FINISH);
    INSERT_ELEMENT(ALARM_SCHEDULER);
#undef INSERT_ELEMENT
  }
  return out << strings[value];
}

bool
TimeLine::AddEvent(Event *ePtr) {
  if (!ePtr || ePtr->GetEventTime() < 0) {
    return false;
  }
  double tm = ePtr->GetEventTime();
  vector<Event *>::iterator it;
  for (it = m_timeline.begin(); it != m_timeline.end(); it++) {
    if (tm < (*it)->GetEventTime()) {
      m_timeline.insert(it, ePtr);
      break;
    }
  }
  if (it == m_timeline.end()) {
    m_timeline.push_back(ePtr);
  }
//    cout << string(FLOAT_TIME_WIDTH+2, ' ')
//         << "[TimeLine::AddEvent] add event "
//         << ePtr->GetEventType() << " at "
//         << ePtr->GetEventTime() << endl;
  return true;
}

Event *
TimeLine::PeekNext() {
  if (!m_timeline.empty()) {
    return m_timeline.front();
  } else {
    cout << "[TimeLine::PeekNext] Error: timeline empty" << endl;
    return NULL;
  }
}

Event *
TimeLine::PopNext() {
  if (!m_timeline.empty()) {
    Event *ePtr2pop = m_timeline.front();
    m_timeline.erase(m_timeline.begin());
    return ePtr2pop;
  } else {
    cout << "[TimeLine::PopNext] Error: timeline empty" << endl;
    return NULL;
  }
}

bool
TimeLine::isEmpty() {
  return m_timeline.empty();
}

Simulator::Simulator() {

  m_currentTime = 0;

  m_schedulerPtr = NULL;
  m_trafficPtr = NULL;
}

bool
Simulator::InstallScheduler(string schedulerName) {
  // can only be installed once
  if (m_schedulerPtr) return false;
  if (0 == schedulerName.compare("sunflow")) {
    m_schedulerPtr = new SchedulerSunflow();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  } else if (0 == schedulerName.compare("Solstice")) {
    m_schedulerPtr = new SchedulerSolstice();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  } else if (0 == schedulerName.compare("BvN")) {
    m_schedulerPtr = new SchedulerBvN();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  } else if (0 == schedulerName.compare("TMS")) {
    m_schedulerPtr = new SchedulerTMS();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  } else if (0 == schedulerName.compare("Edmond")) {
    m_schedulerPtr = new SchedulerEdmond();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  } else if (0 == schedulerName.compare("varysPaper")) {
    m_schedulerPtr = new SchedulerVarysPaper();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  } else if (0 == schedulerName.compare("varysImpl")) {
    m_schedulerPtr = new SchedulerVarysImpl();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  } else if (0 == schedulerName.compare("varysImpl201606")) {
    // ENABLE_PERTURB_IN_PLAY = false;
    // EQUAL_FLOW_TO_SAME_REDUCER = true;
    m_schedulerPtr = new SchedulerVarysImpl201606();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  } else if (0 == schedulerName.compare("aaloImpl")) {
    m_schedulerPtr = new SchedulerAaloImpl();
    m_schedulerPtr->InstallSimulator(this);
    return true;
  }

  return false;
}

bool
Simulator::InstallTrafficGen(string trafficProducerName) {
  // can only be installed once
  if (m_trafficPtr) return false;
  if (0 == trafficProducerName.compare("fb")) {
    m_trafficPtr = new TGTraceFB();
    m_trafficPtr->InstallSimulator(this);
    return true;
  } else if (0 == trafficProducerName.compare("utilizationONLY")) {
    m_trafficPtr = new TGFBUtilizationOnly();
    m_trafficPtr->InstallSimulator(this);
    return true;
  } else if (0 == trafficProducerName.compare("analyzeONLY")) {
    m_trafficPtr = new TGFBAnalyzeOnly();
    m_trafficPtr->InstallSimulator(this);
    return true;
  } else if (0 == trafficProducerName.compare("fb1by1")) {
    m_trafficPtr = new TGFBOnebyOne();
    m_trafficPtr->InstallSimulator(this);
    // comp time stat only available under "1by1".
    LOG_COMP_STAT = true;
    // Should be false unless hacking to obtain quick results.
    SPEEDUP_1BY1_TRAFFIC_IN_SCHEDULER = false;
    return true;
  }

  return false;
}

Simulator::~Simulator() {
  m_timeline.clear();
  delete m_schedulerPtr;
  delete m_trafficPtr;
}

void
Simulator::UpdateSchedulerAlarm(Event *ep) {

  if (m_timeline.empty()) {
    m_timeline.push_back(ep);
    return;
  }
  // from the event next to the current
  for (vector<Event *>::iterator tlIt = m_timeline.begin() + 1;
       tlIt != m_timeline.end(); tlIt++) {
    if ((*tlIt)->GetEventType() == ALARM_SCHEDULER) {
      delete *tlIt;
      m_timeline.erase(tlIt);
      break;
      // there should be no more than 1 scheduler alarm
    }
  }

  // insert the new event
  AddEvent(ep);
}

void
Simulator::UpdateTrafficAlarm(Event *ep) {
  if (m_timeline.empty()) {
    m_timeline.push_back(ep);
    return;
  }
  // from the event next to the current
  for (vector<Event *>::iterator tlIt = m_timeline.begin() + 1;
       tlIt != m_timeline.end(); tlIt++) {
    if ((*tlIt)->GetEventType() == ALARM_TRAFFIC) {
      delete *tlIt;
      m_timeline.erase(tlIt);
      break;
      // there should be no more than 1 scheduler alarm
    }
  }

  // insert the new event
  AddEvent(ep);
}

void
Simulator::Run() {
  if (!m_schedulerPtr) {
    cout << "[Simulator::Run] Error: "
         << "Please install scheduler first!" << endl;
    return;
  }
  if (!m_trafficPtr) {
    cout << "[Simulator::Run] Error: "
         << "Please install traffic generator first!" << endl;
    return;
  }

  m_trafficPtr->NotifySimStart();

  while (!m_timeline.empty()) {
    // Step 1: peek for time of next event
    Event *nextEvent = PeekNext();
    double nextEventTime = nextEvent->GetEventTime();

    if (nextEventTime < m_currentTime) {
      cout << "Error: event time is earlier" <<
           " than current time!" << endl;
    }

    // Step 2: adjust clock
    m_currentTime = nextEventTime;

    // Step 4: execute next event
    switch (nextEvent->GetEventType()) {
      case ALARM_TRAFFIC:DoAlarmTrafficGen();
        break;
      case MSG_TRAFFIC_FINISH:DoNotifyTrafficFinish();
        break;
      case ALARM_SCHEDULER:DoAlarmScheduler();
        break;
      case MSG_ADD_FLOWS:DoNotifyAddFlows();
        break;
      case MSG_ADD_COFLOWS:DoNotifyAddCoflows();
        break;
      default:break;
    }

    // Step 5: rm event executed
    Event *e2p = PopNext();
    if (e2p != nextEvent) {
      cout << "[Simulator::Run] error: "
           << " e2p != currentEvent " << endl;
    }
    delete nextEvent;
  }
  // end of simulation.
  m_trafficPtr->NotifySimEnd();
  m_schedulerPtr->NotifySimEnd();
}

void
Simulator::DoAlarmTrafficGen() {
  m_trafficPtr->TrafficGenAlarmPortal(m_currentTime);
}

void
Simulator::DoNotifyTrafficFinish() {
  //m_trafficPtr->NotifyTrafficReq();
  MsgEventTrafficFinish
      *trafficFinishEvent = (MsgEventTrafficFinish *) PeekNext();
  if (trafficFinishEvent->GetEventType() != MSG_TRAFFIC_FINISH) {
    cout << "[Simulator::DoNotifyTrafficFinish] error: "
         << " the event type is not MSG_TRAFFIC_FINISH!" << endl;
    return;
  }

  m_trafficPtr->NotifyTrafficFinish(m_currentTime,
                                    trafficFinishEvent->m_cfpVp,
                                    trafficFinishEvent->m_fpVp);
}

void
Simulator::DoAlarmScheduler() {
  m_schedulerPtr->SchedulerAlarmPortal(m_currentTime);
}

void
Simulator::DoNotifyAddFlows() {
  MsgEventAddFlows *addFlowsEvent = (MsgEventAddFlows *) PeekNext();
  if (addFlowsEvent->GetEventType() != MSG_ADD_FLOWS) {
    cout << "[Simulator::DoNotifyAddFlows] error: "
         << " the event type is not MSG_ADD_FLOWS!" << endl;
    return;
  }
  m_schedulerPtr->NotifyAddFlows(m_currentTime);
}

void
Simulator::DoNotifyAddCoflows() {
  //m_schedulerPtr->NotifyScheduleEnd();
  MsgEventAddCoflows *addCoflowsEvent = (MsgEventAddCoflows *) PeekNext();
  if (addCoflowsEvent->GetEventType() != MSG_ADD_COFLOWS) {
    cout << "[Simulator::DoNotifyAddCoflows] error: "
         << " the event type is not MSG_ADD_COFLOWS!" << endl;
    return;
  }
  m_schedulerPtr->NotifyAddCoflows(m_currentTime, addCoflowsEvent->m_cfpVp);

}

bool
SchedulerTimeLine::RemoveMultipleEvent(EventType tp) {
  if (m_timeline.size() <= 0) {
    return true;
  }
  if (tp != APPLY_CIRCUIT) {
    // sanity check the legal singular event type
    cout << string(FLOAT_TIME_WIDTH + 2, ' ')
         << "[Simulator::RemoveSingularEvent] ERROR: try to remove " << tp
         << "but it can NOT be multiple!" << endl;
    return false;
  }
  // skip the event we are handling this type of event
  for (vector<Event *>::iterator
           it = m_timeline.begin() + 1; /*remove any NON-CURRENT*/
       it != m_timeline.end();) {
    if ((*it)->GetEventType() == tp) {
      delete *it;
      it = m_timeline.erase(it);
    } else {
      it++;
    }
  }
  return true;
}

bool
SchedulerTimeLine::RemoveSingularEvent(EventType tp) {
  if (m_timeline.size() <= 0) {
    return true;
  }
  if (tp != FLOW_FINISH
      && tp != APPLY_CIRCUIT
      && tp != ORDER_CIRCUIT
      && tp != APPLY_NEW_SCHEDULE
      && tp != RESCHEDULE
      && tp != SCHEDULE_END) {
    // sanity check the legal singular event type
    cout << string(FLOAT_TIME_WIDTH + 2, ' ')
         << "[Simulator::RemoveSingularEvent] ERROR: try to remove " << tp
         << "but it is NOT singular!" << endl;
    return false;
  }
  vector<Event *>::iterator it;
  if ((*m_timeline.begin())->GetEventType() == tp) {
    // we are now handling this type of event
    // there should not be other event with type *tp*
    // in the event vector
    return true;
  }
  for (it = m_timeline.begin() + 1; /*remove any NON-CURRENT*/
       it != m_timeline.end();) {
    if ((*it)->GetEventType() == tp) {
      // important: delete event pointer
      delete *it;
      it = m_timeline.erase(it);
      // there should be no more than one event with type
      // tp in the event vector
      break;
    } else {
      it++;
    }
  }
  return true;
}

void
TimeLine::Print(void) {
  int eventCount = 0;
  for (vector<Event *>::iterator it = m_timeline.begin();
       it != m_timeline.end(); it++) {
    cout << string(FLOAT_TIME_WIDTH + 2, ' ')
         << "[TimeLine::Print] "
         << "TimeLine(" << eventCount++ << ") " << " to do "
         << (*it)->GetEventType() << " at " << (*it)->GetEventTime() << "s"
         << endl;
  }
}