//
//  scheduler.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 9/21/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//

#include <algorithm>
#include <iomanip>
#include <cfloat> 
#include <sys/time.h>
#include <string.h>

#include "scheduler.h"
#include "trafficGen.h"
#include "events.h"
#include "global.h"
#include "util.h"
#include "coflow.h"

#define MWM_RANGE 100000000 //2147483647 = 2,147,483,647

using namespace std;

///////////////////////////////////////////////////////
////////////// Code for base class Scheduler
///////////////////////////////////////////////////////

Scheduler::Scheduler() {
  m_simPtr = NULL;
  m_currentTime = 0;
  m_myTimeLine = new SchedulerTimeLine();
  m_coflowPtrVector = vector<Coflow *>();

  m_nextElecRate = map<int, long>();
  m_nextOptcRate = map<int, long>();

  m_txAuditFile.open(NET_TX_AUDIT_FILE_NAME);
  if (!m_txAuditFile.is_open()) {
    cout << "Error: unable to open network audit file "
         << NET_TX_AUDIT_FILE_NAME << endl;
    cout << "Now terminate the program" << endl;
    exit(-1);
  }
  m_txAuditFile << "time num_coflow num_flow optc_rate(Mbs) elec_rate(Mbs)"
                << endl;

  m_portAuditFile.open(PORT_AUDIT_FILE_NAME);
  if (!m_portAuditFile.is_open()) {
    cout << "Error: unable to open port audit file "
         << PORT_AUDIT_FILE_NAME << endl;
    cout << "Now terminate the program" << endl;
    exit(-1);
  }
  m_port_audit_title_line_last_gone = false;
  m_port_audit_title_line_bottleneck = false;
  m_port_audit_title_line_all = false;
}

Scheduler::~Scheduler() {
  delete m_myTimeLine;

  if (m_txAuditFile.is_open()) {
    m_txAuditFile.close();
  }

  if (m_portAuditFile.is_open()) {
    m_portAuditFile.close();
  }
}

void
Scheduler::UpdateAlarm() {
  // update alarm for scheduler on simulator
  if (!m_myTimeLine->isEmpty()) {
    Event *nextEvent = m_myTimeLine->PeekNext();
    double nextTime = nextEvent->GetEventTime();
    Event *schedulerAlarm = new Event(ALARM_SCHEDULER, nextTime);
    m_simPtr->UpdateSchedulerAlarm(schedulerAlarm);
  }
}

void
Scheduler::NotifySimEnd() {
  return;
  cout << "[Scheduler::NotifySimEnd()] is called." << endl;
  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {

      if ((*fpIt)->GetBitsLeft() <= 0) {
        // such flow has finished
        continue;
      }

      cout << "[Scheduler::NotifySimEnd()] flow [" << (*fpIt)->GetFlowId()
           << "] "
           << "(" << (*fpIt)->GetSrc() << "=>" << (*fpIt)->GetDest() << ") "
           << (*fpIt)->GetSizeInBit() << " bytes "
           << (*fpIt)->GetBitsLeft() << " bytes left "
           << (*fpIt)->GetElecRate() << " bps" << endl;

    }
  }
}

bool
Scheduler::Transmit(double startTime,
                    double endTime,
                    bool basic,
                    bool local,
                    bool salvage) {

  CircuitAuditIfNeeded(startTime, endTime);

  bool hasCoflowFinish = false;
  bool hasFlowFinish = false;
  bool hasCoflowTmpFinish = false;

  vector<Coflow *> finished_coflows;
  vector<Flow *> finished_flows;

  m_validate_last_tx_src_bits.clear();
  m_validate_last_tx_dst_bits.clear();

  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end();) {

    if ((*cfIt)->IsFlowsAddedComplete()) {
      // coflow not completed yet
      // but the flows added so far have finished
      cfIt++;
      continue;
    }

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {

      if ((*fpIt)->GetBitsLeft() <= 0) {
        // such flow has finished
        continue;
      }

      // tx rate verification debug
      long validate_tx_this_flow_bits = (*fpIt)->GetBitsLeft();

      // ********* begin tx ****************
      if (basic) {
        (*fpIt)->Transmit(startTime, endTime);
      }
      if (local) {
        (*fpIt)->TxLocal();
      }
      if (salvage) {
        (*fpIt)->TxSalvage();
      }
      // ********* end tx ********************

      // tx rate verification debug
      validate_tx_this_flow_bits -= (*fpIt)->GetBitsLeft();
      MapWithInc(m_validate_last_tx_src_bits,
                 (*fpIt)->GetSrc(),
                 validate_tx_this_flow_bits);
      MapWithInc(m_validate_last_tx_dst_bits,
                 (*fpIt)->GetDest(),
                 validate_tx_this_flow_bits);

      if ((*fpIt)->GetBitsLeft() == 0) {
        hasFlowFinish = true;
        (*cfIt)->NumFlowFinishInc();
        (*fpIt)->SetEndTime(endTime);
        finished_flows.push_back(*fpIt);
      }

      // update coflow account on bytes sent.
      (*cfIt)->AddTxBit(validate_tx_this_flow_bits);
    }

    // debug for coflow progress
    if (DEBUG_LEVEL >= 3 && hasFlowFinish) {
      cout << fixed << setw(FLOAT_TIME_WIDTH) << endTime << "s ";
      cout << (*cfIt)->toString() << endl;
    }

    if ((*cfIt)->IsComplete()) {

      //cout << string(FLOAT_TIME_WIDTH+2, ' ')
      cout << fixed << setw(FLOAT_TIME_WIDTH) << endTime << "s "
           << "[Scheduler::Transmit] coflow finish! # "
           << (*cfIt)->GetJobId() << endl;

      finished_coflows.push_back(*cfIt);
      (*cfIt)->SetEndTime(endTime);
      cfIt = m_coflowPtrVector.erase(cfIt);
      hasCoflowFinish = true;

    } else {

      if ((*cfIt)->IsFlowsAddedComplete()) {
        hasCoflowTmpFinish = true;
      }
      // jump to next coflow
      cfIt++;
    }
  }

  ScheduleToNotifyTrafficFinish(endTime, finished_coflows, finished_flows);

  if (hasCoflowFinish || hasCoflowTmpFinish) {
    CoflowFinishCallBack(endTime);
  } else if (hasFlowFinish) {
    FlowFinishCallBack(endTime);
  }

  if (hasCoflowFinish || hasCoflowTmpFinish || hasFlowFinish) {
    Scheduler::UpdateFlowFinishEvent(endTime);
  }

  if (hasCoflowFinish || hasCoflowTmpFinish || hasFlowFinish) {
    return true;
  }
  return false;
}

void
Scheduler::ScheduleToNotifyTrafficFinish(double end_time,
                                         vector<Coflow *> &coflows_done,
                                         vector<Flow *> &flows_done) {
  if (coflows_done.empty() && flows_done.empty()) {
    return;
  }

  if (coflows_done.empty()
      && SPEEDUP_1BY1_TRAFFIC_IN_SCHEDULER
      && m_coflowPtrVector.size() == 1) {
    // no FCT for hacking mode.
    return;
  }
  //notify traffic generator of coflow / flow finish
  vector<Coflow *> *finishedCf = new vector<Coflow *>(coflows_done);
  vector<Flow *> *finishedF = new vector<Flow *>(flows_done);
  MsgEventTrafficFinish *msgEventPtr =
      new MsgEventTrafficFinish(end_time, finishedCf, finishedF);
  m_simPtr->AddEvent(msgEventPtr);

  if (!coflows_done.empty()) {
    WriteCircuitAuditIfNeeded(end_time, coflows_done);
  }
}

//returns negative if all flows has finished
//return DBL_MAX if all flows are waiting indefinitely
double
Scheduler::CalcTime2FirstFlowEnd() {
  double time2FirstFinish = DBL_MAX;
  bool hasUnfinishedFlow = false;
  bool finishTimeValid = false;
  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {

    if ((*cfIt)->IsFlowsAddedComplete()) {
      //flows added in such coflow have all completed
      continue;
    }

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {

      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }

      hasUnfinishedFlow = true;

      // calc the min finishing time
      double flowCompleteTime = DBL_MAX;
      if ((*fpIt)->isThruOptic() && (*fpIt)->GetOptcRate() > 0) {
        flowCompleteTime = SecureFinishTime((*fpIt)->GetBitsLeft(),
                                            (*fpIt)->GetOptcRate());
      } else if (!(*fpIt)->isThruOptic() && (*fpIt)->GetElecRate() > 0) {
        flowCompleteTime = SecureFinishTime((*fpIt)->GetBitsLeft(),
                                            (*fpIt)->GetElecRate());
      }
      if (time2FirstFinish > flowCompleteTime) {
        finishTimeValid = true;
        time2FirstFinish = flowCompleteTime;
      }
    }
  }
  if (hasUnfinishedFlow) {
    if (finishTimeValid) {
      return time2FirstFinish;
    } else {
      //all flows are waiting indefinitely
      return DBL_MAX;
    }
  } else {
    // all flows are finished
    return -DBL_MAX;
  }
}

void
Scheduler::UpdateFlowFinishEvent(double baseTime) {

  double time2FirstFinish = CalcTime2FirstFlowEnd();

  if (time2FirstFinish == DBL_MAX) {
    // all flows are waiting indefinatly
    m_myTimeLine->RemoveSingularEvent(FLOW_FINISH);

  } else if (time2FirstFinish == -DBL_MAX) {
    // all flows are done

  } else {
    // valid finishing time
    m_myTimeLine->RemoveSingularEvent(FLOW_FINISH);
    double firstFinishTime = baseTime + time2FirstFinish;
    Event *flowFinishEventPtr = new Event(FLOW_FINISH, firstFinishTime);
    m_myTimeLine->AddEvent(flowFinishEventPtr);
  }
}

void
Scheduler::UpdateRescheduleEvent(double reScheduleTime) {
  m_myTimeLine->RemoveSingularEvent(RESCHEDULE);
  Event *rescheduleEventPtr = new Event(RESCHEDULE, reScheduleTime);
  m_myTimeLine->AddEvent(rescheduleEventPtr);
  /*
  cout << string(2+FLOAT_TIME_WIDTH, ' ')
  << "[Scheduler::UpdateRescheduleEvent] try to reschedule at "
  << reScheduleTime << "s" <<endl;
   */
}

void
Scheduler::NotifyAddFlows(double alarmTime) {
  //FlowArrive(alarmTime);
  EventFlowArrive *msgEp = new EventFlowArrive(alarmTime);
  m_myTimeLine->AddEvent(msgEp);
  UpdateAlarm();
}

void
Scheduler::NotifyAddCoflows(double alarmTime, vector<Coflow *> *cfVecPtr) {
  //CoflowArrive(alarmTime,cfVecPtr);
  EventCoflowArrive *msgEp = new EventCoflowArrive(alarmTime, cfVecPtr);
  m_myTimeLine->AddEvent(msgEp);
  UpdateAlarm();
}

double
Scheduler::SecureFinishTime(long bits, long rate) {
  if (rate == 0) {
    return DBL_MAX;
  }
  double timeLen = (double) bits / (double) rate;

  // more complicated impl below *********
  long bitsleft = 1;
  int delta = 0;
  while (bitsleft > 0) {
    timeLen = (double) (delta + bits) / (double) rate;
    bitsleft = bits - rate * timeLen;
    delta++;
  }
  // more complicated impl above *****

  //if (timeLen < 0.0000000001) return 0.0000000001;
  return timeLen;
}

void
Scheduler::Print(void) {
  return;
  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {
    if (cfIt == m_coflowPtrVector.begin()) {
      cout << fixed << setw(FLOAT_TIME_WIDTH)
           << m_currentTime << "s ";
    } else {
      cout << string(FLOAT_TIME_WIDTH + 2, ' ');
    }

    cout << "[Scheduler::Print] "
         << "Coflow ID " << (*cfIt)->GetCoflowId() << endl;
    (*cfIt)->Print();
  }
}

// copy flow rate from m_nextElecRate & m_nextOptcRate
// and reflect the rate on flow record.
void
Scheduler::SetFlowRate() {
  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {
    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();
    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      //set flow rate
      int flowId = (*fpIt)->GetFlowId();
      long elecBps = MapWithDef(m_nextElecRate, flowId, (long) 0);
      long optcBps = MapWithDef(m_nextOptcRate, flowId, (long) 0);
      (*fpIt)->SetRate(elecBps, optcBps);
    }
  }
}

void
Scheduler::CalAlphaAndSortCoflowsInPlace(vector<Coflow *> &coflows) {
  for (vector<Coflow *>::iterator it = coflows.begin();
       it != coflows.end(); it++) {
    (*it)->CalcAlpha();
  }
  std::stable_sort(coflows.begin(), coflows.end(), coflowCompAlpha);
}

void
Scheduler::sortCoflowOnStaticAlphaInPlace(vector<Coflow *> &coflows) {
  std::stable_sort(coflows.begin(), coflows.end(), coflowCompAlpha);
}

void
Scheduler::SortCoflowsInPlaceArrival(vector<Coflow *> &coflows) {
  std::stable_sort(coflows.begin(), coflows.end(), coflowCompArrival);
}

bool
Scheduler::ValidateLastTxMeetConstraints(long port_bound_bits) {
  for (map<int, long>::const_iterator
           src_kv_pair = m_validate_last_tx_src_bits.begin();
       src_kv_pair != m_validate_last_tx_src_bits.end();
       src_kv_pair++) {
    if (src_kv_pair->second > port_bound_bits) {
      cout << "Error in validating TX constraints!!! " << endl
           << " src " << src_kv_pair->first << " flows over bound "
           << port_bound_bits << endl;
      return false;
    }
  }

  for (map<int, long>::const_iterator
           dst_kv_pair = m_validate_last_tx_dst_bits.begin();
       dst_kv_pair != m_validate_last_tx_dst_bits.end();
       dst_kv_pair++) {
    if (dst_kv_pair->second > port_bound_bits) {
      cout << "Error in validating TX constraints!!! " << endl
           << " dst " << dst_kv_pair->first << "  flows over bound "
           << port_bound_bits << endl;
      return false;
    }
  }
  // cout << "TX budget is valid." << endl;
  return true;
}
