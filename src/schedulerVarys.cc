//
//  schedulerVarys.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 2/5/15.
//  Copyright (c) 2015 Xin Sunny Huang. All rights reserved.
//

#include <algorithm>
#include <iomanip>
#include <sys/time.h>

#include "scheduler.h"
#include "events.h"
#include "global.h"
#include "util.h"
#include "coflow.h" 

///////////////////////////////////////////////////////
////////////// Code for Varys
///////////////////////////////////////////////////////

SchedulerVarys::SchedulerVarys() : Scheduler() {

}

SchedulerVarys::~SchedulerVarys() {
}

void
SchedulerVarys::SchedulerAlarmPortal(double alarmTime) {

  while (!m_myTimeLine->isEmpty()) {
    Event *currentEvent = m_myTimeLine->PeekNext();
    double currentEventTime = currentEvent->GetEventTime();

    if (currentEventTime > alarmTime) {
      break;
    }

//        cout << fixed << setw(FLOAT_TIME_WIDTH) << currentEventTime << "s "
//        << "[SchedulerVarys::SchedulerAlarmPortal] "
//        << " working on event type " << currentEvent->GetEventType() << endl;

    bool has_flow_finished = false;
    if (m_currentTime < alarmTime) {
      has_flow_finished = Transmit(m_currentTime, alarmTime,
                                   true,  // basic
                                   true,  // local - should be useless
                                   FLOW_FINISH
                                       == currentEvent->GetEventType() // salvage
      );
      m_currentTime = alarmTime;
    }

    switch (currentEvent->GetEventType()) {
      case RESCHEDULE:Schedule();
        break;
      case COFLOW_ARRIVE:CoflowArrive();
        // clean any local traffic
        Transmit(m_currentTime,
                 m_currentTime, /*basic*/
                 false, /*local*/
                 true, /*salvage*/
                 false);
        break;
      case FLOW_ARRIVE:FlowArrive();
        break;
      case APPLY_NEW_SCHEDULE:ApplyNewSchedule();
        break;
      case FLOW_FINISH:break;
      default:break;
    }
    Event *e2p = m_myTimeLine->PopNext();
    if (e2p != currentEvent) {
      cout << "[SchedulerVarys::SchedulerAlarmPortal] error: "
           << " e2p != currentEvent " << endl;
    }
    delete currentEvent;
  }

  UpdateAlarm();
}

void
SchedulerVarys::ApplyNewSchedule() {
  // Capture changes in rate allocation.
  // 1. Reflect the rate in m_nextElecRate & m_nextOptcRate to flow record.
  Scheduler::SetFlowRate();
  // 2. Update flow finish event depending on the latest rate && circuit.
  Scheduler::UpdateFlowFinishEvent(m_currentTime);
}

void
SchedulerVarys::CoflowArrive() {
  // unbox coflow vector pointer
  EventCoflowArrive
      *coflowsArriveEvent = (EventCoflowArrive *) m_myTimeLine->PeekNext();
  if (coflowsArriveEvent->GetEventType() != COFLOW_ARRIVE) {
    cout << "[SchedulerVarys::CoflowArrive] error: "
         << " the event type is not COFLOW_ARRIVE!" << endl;
    return;
  }
  AddCoflows(coflowsArriveEvent->m_cfpVp);
  Scheduler::UpdateRescheduleEvent(m_currentTime);
}

void
SchedulerVarys::AddCoflows(vector<Coflow *> *cfsPtr) {
  if (!cfsPtr) return;
  // check all coflows are legal
  for (vector<Coflow *>::iterator cfpIt = cfsPtr->begin();
       cfpIt != cfsPtr->end(); cfpIt++) {
    if (!(*cfpIt)->isRawCoflow()) {
      cout << "[SchedulerVarys::AddCoflows]"
           << "error: is not a raw coflow" << endl;
      continue;
    }
    m_coflowPtrVector.push_back(*cfpIt);
  }
}

void
SchedulerVarys::FlowArrive() {
}

void
SchedulerVarys::AddFlows() {
}

void
SchedulerVarys::CoflowFinishCallBack(double finishTime) {
  Scheduler::UpdateRescheduleEvent(finishTime);
}

void
SchedulerVarys::FlowFinishCallBack(double finishTime) {
}


//
// SchedulerVarysImpl :
//          Simulator implementation as seen on github at some time in 2015
//

void
SchedulerVarysImpl::Schedule() {

  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerVarysImpl::schedule] VarysImpl scheduling START" << endl;

  Print();

  struct timeval start_time;
  gettimeofday(&start_time, NULL);
  /////////////// varys //////////////////////////


  // STEP 1: Initialize next rate for all flows to (0,0)
  m_nextElecRate.clear();
  m_nextOptcRate.clear();

  // STEP 2: Perform varys rate control
  RateControlVarysImpl(m_coflowPtrVector, m_nextElecRate, ELEC_BPS);

  struct timeval end_time;
  gettimeofday(&end_time, NULL);

  double ComputationSeconds = secondPass(end_time, start_time);
  // test of pattern computation delay
  if (ZERO_COMP_TIME) {
    ComputationSeconds = 0.0;
  }

  //debug
  /*
   cout << fixed << setw(FLOAT_TIME_WIDTH)
   << m_currentTime << "s "
   << "[SchedulerVarysImpl::schedule] VarysImpl scheduling DONE in "
   << ComputationSeconds << "s"<< endl;
   */

  m_myTimeLine->RemoveSingularEvent(APPLY_NEW_SCHEDULE);
  double activateTime = m_currentTime + ComputationSeconds;
  Event *applyScheduleEventPtr = new Event(APPLY_NEW_SCHEDULE, activateTime);
  m_myTimeLine->AddEvent(applyScheduleEventPtr);

  //debug
  /*
   cout << fixed << setw(FLOAT_TIME_WIDTH)
   << m_currentTime << "s "
   << "[SchedulerVarysImpl::schedule] new schedule will be activated at "
   << activateTime  << "s" << endl;
   */
}

// perform varys based rate control - a similar version as seen in Github some time in 2015.
// sort coflows (in place).
// recored flow rate to rates (in place).
// coflows: a list of coflows
// rates : map from flow id to allocated rate.
void
SchedulerVarysImpl::RateControlVarysImpl(vector<Coflow *> &coflows,
                                         map<int, long> &rates,
                                         long LINK_RATE_BPS) {

  if (coflows.empty()) {
    return;
  }

  rates.clear();

  // STEP 1: Sort ALL coflows based on different scheduling policies.
  // varys reschedules upon coflow arrival/departure =>
  // coflows are sorted upon arrival/departure.

  // fixed by Sunny : sort Coflows on runtime bottle-neck.
  // The performance is much better than original implementation.
  CalAlphaAndSortCoflowsInPlace(coflows);

  // initialize.
  map<int, long> sBpsFree, rBpsFree;

  for (unsigned int first_unscheduled_cf_idx = 0;
       first_unscheduled_cf_idx < coflows.size();
       first_unscheduled_cf_idx++) {

    Coflow *cf_to_be_schedule = coflows[first_unscheduled_cf_idx];

    // for each coflow
    if (cf_to_be_schedule->IsComplete()
        || 0 == cf_to_be_schedule->GetAlpha()) {
      //such coflow has completed
      // or has zero demand
      continue;
    }

    vector<Flow *> *flowVecPtr = cf_to_be_schedule->GetFlows();

    map<int, long> sBpsUsed, rBpsUsed;

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }
      /*
       // In the original implementation, residual bandwidth is updated
       // after rate-allocation for each coflow.
       // Therefore some of the coflow that has a conflict with the prioritized coflows
       // (i.e. the less prioritized coflow request a link that is fully utilized by)
       // a previous coflow. may receive 0 allocated rate.
       long sBpsMax = MapWithDef(m_varys_sBpsFree, (*fpIt)->GetSrc(), LINK_RATE_BPS);
       long rBpsMax = MapWithDef(m_varys_rBpsFree, (*fpIt)->GetDest(), LINK_RATE_BPS);
       long feasibleBpsMax = sBpsMax < rBpsMax ? sBpsMax : rBpsMax;
       if (flowBitps > feasibleBpsMax) {
       flowBitps = feasibleBpsMax;
       }
       */

      // another proposal - more selfish coflow.
      // this proposal has much better performance.
      long sBps = MapWithDef(sBpsFree, (*fpIt)->GetSrc(), LINK_RATE_BPS);
      long rBps = MapWithDef(rBpsFree, (*fpIt)->GetDest(), LINK_RATE_BPS);
      long minFreeBps = sBps < rBps ? sBps : rBps;
      // Intuition is as follows:
      // Assume the current flow is on the bottleneck port.
      // Therefore the bottleneck link rate is minFreeBps.
      // Hence, according to MADD, allocate rate based on data size.
      // Note that flowBitps might be less than it should be.
      // because alpha is the max estimate of the
      // max( sum(src-demand-of-this-coflow),
      //      sum(dst-demand-of-this-coflow))
      long flowBitps = minFreeBps * ((*fpIt)->GetBitsLeft()
          / (double) cf_to_be_schedule->GetAlpha());


      // update utilization profile.
      if (flowBitps > 0) {
        MapWithDef(rates, (*fpIt)->GetFlowId(), flowBitps);
        MapWithInc(sBpsUsed, (*fpIt)->GetSrc(), flowBitps);
        MapWithInc(rBpsUsed, (*fpIt)->GetDest(), flowBitps);
      }
    } // for each flow

    // Remove capacity from ALL sources and destination for this coflow
    for (map<int, long>::iterator sUsedIt = sBpsUsed.begin();
         sUsedIt != sBpsUsed.end(); sUsedIt++) {
      MapWithDef(sBpsFree, sUsedIt->first, LINK_RATE_BPS);
      sBpsFree[sUsedIt->first] -= sUsedIt->second;
    }
    for (map<int, long>::iterator rUsedIt = rBpsUsed.begin();
         rUsedIt != rBpsUsed.end(); rUsedIt++) {
      MapWithDef(rBpsFree, rUsedIt->first, LINK_RATE_BPS);
      rBpsFree[rUsedIt->first] -= rUsedIt->second;
    }
  } // for each coflow.

  m_sBpsFree_before_workConserv.clear();
  m_sBpsFree_before_workConserv = sBpsFree;
  m_rBpsFree_before_workConserv.clear();
  m_rBpsFree_before_workConserv = rBpsFree;

  // STEP2A: Work conservation as seen in Github.
  RateControlWorkConservationImpl(coflows, rates,
                                  sBpsFree, rBpsFree,
                                  LINK_RATE_BPS);
}

// perform work conservation in the order of coflows, and flows within a coflow.
// given available src/dst port bandwidth resource left
// in sBpsFree & rBpsFree.
// specify link bandwidth capacity as LINK_RATE_BPS.
// record rate in rates.
// NOTE: rates may be modified in place.
void
SchedulerVarysImpl::RateControlWorkConservationImpl(vector<Coflow *> &coflows,
                                                    map<int, long> &rates,
                                                    map<int, long> &sBpsFree,
                                                    map<int, long> &rBpsFree,
                                                    long LINK_RATE_BPS) {

  // Original heuristic: Sort coflows by arrival time and then refill.
  // fixed by Sunny: if not in deadline mode, just maintain the order based on
  // run-time bottle-neck. The performance is much better than original
  // implementation.
  vector<Coflow *> sortedByEDF = coflows;
  // SortCoflowsInPlaceEDFIfNeeded(sortedByEDF);

  for (vector<Coflow *>::iterator cfIt = sortedByEDF.begin();
       cfIt != sortedByEDF.end(); cfIt++) {

    if ((*cfIt)->IsComplete()) {
      //such coflow has completed
      continue;
    }

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {

      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }
      long sBps = MapWithDef(sBpsFree, (*fpIt)->GetSrc(), LINK_RATE_BPS);
      long rBps = MapWithDef(rBpsFree, (*fpIt)->GetDest(), LINK_RATE_BPS);
      long minFreeBps = sBps < rBps ? sBps : rBps;
      if (minFreeBps > 0) {
        MapWithInc(rates, (*fpIt)->GetFlowId(), minFreeBps);
        sBpsFree[(*fpIt)->GetSrc()] -= minFreeBps;
        rBpsFree[(*fpIt)->GetDest()] -= minFreeBps;
      }
    }
  }

  //debug
  if (DEBUG_LEVEL >= 5) {
    map<int, long> rates_debug = rates;
    cout << " demand " << endl;
    for (vector<Coflow *>::const_iterator
             cf_iter = sortedByEDF.begin();
         cf_iter != sortedByEDF.end();
         cf_iter++) {

      Coflow *cf = *cf_iter;
      cout << " coflow id " << cf->GetCoflowId()
           << " alpha " << cf->GetAlpha()
           << " online_alpha " << cf->GetAlphaOnline() << endl;

      vector<Flow *> flows = *(cf->GetFlows());
      for (vector<Flow *>::const_iterator f_iter = flows.begin();
           f_iter != flows.end(); f_iter++) {
        Flow *flow = *f_iter;
        cout << " flow id [" << flow->GetFlowId() << "] "
             << flow->GetSrc() << "->" << flow->GetDest()
             << " demand " << flow->GetBitsLeft()
             << " rate " << MapWithDef(rates_debug, flow->GetFlowId(), (long) 0)
             << endl;
      }
    }
  }
}

//
// SchedulerVarysImpl201606 :
// Simulation implementation as seen on GitHub in June 2016
//

void
SchedulerVarysImpl201606::Schedule() {

  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerVarysImpl201606::schedule] VarysImpl201606 scheduling START"
       << endl;

  Print();

  struct timeval start_time;
  gettimeofday(&start_time, NULL);

  // STEP 1: Initialize next rate for all flows to (0,0)
  m_nextElecRate.clear();
  m_nextOptcRate.clear();

  // STEP 2: Perform varys rate control
  RateControlVarysImpl201606(m_coflowPtrVector, m_nextElecRate, ELEC_BPS);

  struct timeval end_time;
  gettimeofday(&end_time, NULL);

  double ComputationSeconds = secondPass(end_time, start_time);

  if (ZERO_COMP_TIME) {
    ComputationSeconds = 0.0;
  }

  m_myTimeLine->RemoveSingularEvent(APPLY_NEW_SCHEDULE);
  double activateTime = m_currentTime + ComputationSeconds;
  Event *applyScheduleEventPtr = new Event(APPLY_NEW_SCHEDULE, activateTime);
  m_myTimeLine->AddEvent(applyScheduleEventPtr);
}

void
SchedulerVarysImpl201606::RateControlVarysImpl201606(vector<Coflow *> &coflows,
                                                     map<int, long> &rates,
                                                     long LINK_RATE_BPS) {

  // Mosharaf's implementation is based on synthetic environment:
  // no perturb when generating traffic.

  if (coflows.empty()) {
    return;
  }

  rates.clear();

  // Mosharaf's origninal : sort on static alpha
  sortCoflowOnStaticAlphaInPlace(coflows);

  // initialize.
  map<int, long> sBpsFree, rBpsFree;

  vector<Coflow *> dead_coflows;

  // steps 1: MADD
  for (unsigned int first_unscheduled_cf_idx = 0;
       first_unscheduled_cf_idx < coflows.size();
       first_unscheduled_cf_idx++) {

    Coflow *cf_to_be_schedule = coflows[first_unscheduled_cf_idx];

    // for each coflow
    if (cf_to_be_schedule->IsComplete()) {
      //such coflow has completed
      // or has zero demand
      continue;
    }

    double expected_finish_time_for_cf_to_be_schedule
        = cf_to_be_schedule->CalcAlphaOnline(sBpsFree, rBpsFree, LINK_RATE_BPS);

    if (0 < CF_DEAD_ALPHA_SIGN * expected_finish_time_for_cf_to_be_schedule) {
      dead_coflows.push_back(cf_to_be_schedule);
      continue;//with next coflow
    }

    // valid expected finished time.
    vector<Flow *> *flowVecPtr = cf_to_be_schedule->GetFlows();

    map<int, long> sBpsUsed, rBpsUsed;

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }

      long flowBitps = (*fpIt)->GetBitsLeft()
          / expected_finish_time_for_cf_to_be_schedule;

      long src_cap = MapWithDef(sBpsFree, (*fpIt)->GetSrc(), LINK_RATE_BPS);
      long dst_cap = MapWithDef(rBpsFree, (*fpIt)->GetDest(), LINK_RATE_BPS);

      if (flowBitps > src_cap || flowBitps > dst_cap) {
        cout << "shouldn't be here : what is wrong...?";
        flowBitps = min(src_cap, dst_cap);
      }

      // update utilization profile.
      if (flowBitps > 0) {
        MapWithDef(rates, (*fpIt)->GetFlowId(), flowBitps);
        MapWithInc(sBpsUsed, (*fpIt)->GetSrc(), flowBitps);
        MapWithInc(rBpsUsed, (*fpIt)->GetDest(), flowBitps);
      }
    } // for each flow

    // Remove capacity from ALL sources and destination for this coflow
    for (map<int, long>::iterator sUsedIt = sBpsUsed.begin();
         sUsedIt != sBpsUsed.end(); sUsedIt++) {
      MapWithDef(sBpsFree, sUsedIt->first, LINK_RATE_BPS);
      sBpsFree[sUsedIt->first] -= sUsedIt->second;
    }
    for (map<int, long>::iterator rUsedIt = rBpsUsed.begin();
         rUsedIt != rBpsUsed.end(); rUsedIt++) {
      MapWithDef(rBpsFree, rUsedIt->first, LINK_RATE_BPS);
      rBpsFree[rUsedIt->first] -= rUsedIt->second;
    }
  } // for each coflow.

  // step 2: fair-share among dead coflows which are skipped in round 1 of MADD
  // map from src/dst rack id to number of flows behind the rack
  map<int, int> src_to_num_flow, dst_to_num_flow;
  for (Coflow *one_dead_coflow : dead_coflows) {
    if (one_dead_coflow->IsComplete()) {
      //such coflow has completed
      // or has zero demand
      continue;
    }
    vector<Flow *> *flowVecPtr = one_dead_coflow->GetFlows();
    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }
      MapWithInc(src_to_num_flow, (*fpIt)->GetSrc(), 1);
      MapWithInc(dst_to_num_flow, (*fpIt)->GetDest(), 1);
    }
  } // for each skipped coflow.

  map<int, long> sBpsUsed_fairshare_fix, rBpsUsed_fairshare_fix;
  for (Coflow *one_dead_coflow : dead_coflows) {
    if (one_dead_coflow->IsComplete()) {
      //such coflow has completed
      // or has zero demand
      continue;
    }
    vector<Flow *> *flowVecPtr = one_dead_coflow->GetFlows();
    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }

      long src_bps = MapWithDef(sBpsFree, (*fpIt)->GetSrc(), LINK_RATE_BPS)
          / MapWithDef(src_to_num_flow, (*fpIt)->GetSrc(), 1);
      long dst_bps = MapWithDef(rBpsFree, (*fpIt)->GetDest(), LINK_RATE_BPS)
          / MapWithDef(dst_to_num_flow, (*fpIt)->GetDest(), 1);;

      long flowBitps = min(src_bps, dst_bps);
      // update utilization profile.
      if (flowBitps > 0) {
        MapWithDef(rates, (*fpIt)->GetFlowId(), flowBitps);
        MapWithInc(sBpsUsed_fairshare_fix, (*fpIt)->GetSrc(), flowBitps);
        MapWithInc(rBpsUsed_fairshare_fix, (*fpIt)->GetDest(), flowBitps);
      }
    }
  } // for each skipped coflow. // end of fair share fix.
  // Remove capacity used in the fair share fix.
  for (map<int, long>::iterator sUsedIt = sBpsUsed_fairshare_fix.begin();
       sUsedIt != sBpsUsed_fairshare_fix.end(); sUsedIt++) {
    MapWithDef(sBpsFree, sUsedIt->first, LINK_RATE_BPS);
    sBpsFree[sUsedIt->first] -= sUsedIt->second;
  }
  for (map<int, long>::iterator rUsedIt = rBpsUsed_fairshare_fix.begin();
       rUsedIt != rBpsUsed_fairshare_fix.end(); rUsedIt++) {
    MapWithDef(rBpsFree, rUsedIt->first, LINK_RATE_BPS);
    rBpsFree[rUsedIt->first] -= rUsedIt->second;
  }

  // step 3: Work conservation
  // Heuristic: Sort coflows by arrival time and then refill
  vector<Coflow *> sortedByArrival = coflows;
  SortCoflowsInPlaceArrival(sortedByArrival);

  for (vector<Coflow *>::iterator
           cfIt = sortedByArrival.begin();
       cfIt != sortedByArrival.end();
       cfIt++) {

    if ((*cfIt)->IsComplete()) {
      //such coflow has completed
      continue;
    }

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {

      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }
      long sBps = MapWithDef(sBpsFree, (*fpIt)->GetSrc(), LINK_RATE_BPS);
      long rBps = MapWithDef(rBpsFree, (*fpIt)->GetDest(), LINK_RATE_BPS);
      long minFreeBps = sBps < rBps ? sBps : rBps;
      if (minFreeBps > 0) {
        MapWithInc(rates, (*fpIt)->GetFlowId(), minFreeBps);
        sBpsFree[(*fpIt)->GetSrc()] -= minFreeBps;
        rBpsFree[(*fpIt)->GetDest()] -= minFreeBps;
      }
    }
  }
}

void
SchedulerVarysImpl201606::CheckLinkRateConstraints(const vector<Coflow *> &coflows,
                                                   const map<int, long> &rates,
                                                   long LINK_RATE_BPS) {
  // validate constraints are met.
  map<int, long> varify_src_used, varify_dst_used;
  map<int, long> rates_for_verification = rates;
  for (Coflow *coflow : coflows) {
    for (Flow *flow : *(coflow->GetFlows())) {
      long flow_rate = MapWithDef(rates_for_verification,
                                  flow->GetFlowId(), (long) 0);
      MapWithInc(varify_src_used, flow->GetSrc(), flow_rate);
      MapWithInc(varify_dst_used, flow->GetDest(), flow_rate);
    } // for each flow
  } // for each coflow

  for (auto src_used_pair : varify_src_used) {
    if (src_used_pair.second > LINK_RATE_BPS) {
      cerr << " SRC " << src_used_pair.first
           << " used " << src_used_pair.second
           << " > " << LINK_RATE_BPS;
    }
  }

  for (auto dst_used_pair : varify_dst_used) {
    if (dst_used_pair.second > LINK_RATE_BPS) {
      cerr << " DST " << dst_used_pair.first
           << " used " << dst_used_pair.second
           << " > " << LINK_RATE_BPS;
    }
  }
}

//
//   SchedulerVarys * Varys coflow scheduler As described on paper
//   ERROR : SchedulerVarys just does NOT work.
//

void
SchedulerVarysPaper::Schedule() {

  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerVarysPaper::schedule] VarysPaper scheduling START" << endl;

  Print();

  struct timeval start_time;
  gettimeofday(&start_time, NULL);

  // STEP 1: Initialize next rate for all flows to (0,0)
  m_nextElecRate.clear();
  m_nextOptcRate.clear();

  // STEP 2: Perform varys rate control
  RateControlVarysPaper(m_coflowPtrVector, m_nextElecRate, ELEC_BPS);

  struct timeval end_time;
  gettimeofday(&end_time, NULL);

  double ComputationSeconds = secondPass(end_time, start_time);
  // test of pattern computation delay
  if (ZERO_COMP_TIME) {
    ComputationSeconds = 0.0;
  }

  m_myTimeLine->RemoveSingularEvent(APPLY_NEW_SCHEDULE);
  double activateTime = m_currentTime + ComputationSeconds;
  Event *applyScheduleEventPtr = new Event(APPLY_NEW_SCHEDULE, activateTime);
  m_myTimeLine->AddEvent(applyScheduleEventPtr);

}

// perform varys based rate control - as described in the paper.
// sort coflows (in place).
// recored flow rate to rates (in place).
// coflows: a list of coflows
// rates : map from flow id to allocated rate.
void
SchedulerVarysPaper::RateControlVarysPaper(vector<Coflow *> &coflows,
                                           map<int, long> &rates,
                                           long LINK_RATE_BPS) {

  if (coflows.empty()) {
    return;
  }

  rates.clear();

  // STEP 1: Sort ALL coflows based on different scheduling policies.
  // varys reschedules upon coflow arrival/departure =>
  // coflows are sorted upon arrival/departure.
  CalAlphaAndSortCoflowsInPlace(coflows);

  map<int, long> sBpsFree;
  map<int, long> rBpsFree;

  // MADD herustic (Minimum-Allocation-for-Desired-Duration)
  for (unsigned int first_unscheduled_cf_idx = 0;
       first_unscheduled_cf_idx < coflows.size();
       first_unscheduled_cf_idx++) {

    Coflow *cf_to_be_schedule = coflows[first_unscheduled_cf_idx];

    // for each coflow
    if (cf_to_be_schedule->IsComplete()
        || 0 == cf_to_be_schedule->GetAlpha()) {
      //such coflow has completed
      // or has zero demand
      continue;
    }

    cf_to_be_schedule->CalcAlphaOnline(sBpsFree, rBpsFree, LINK_RATE_BPS);
    if (0 < CF_DEAD_ALPHA_SIGN * cf_to_be_schedule->GetAlphaOnline()) {
      // this coflow can't not be scheduled
      //  because all the available bandwidth have been used up.
      continue;
    }

    vector<Flow *> *flowVecPtr = cf_to_be_schedule->GetFlows();

    map<int, long> sBpsUsed;
    map<int, long> rBpsUsed;

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }

      // MADD proposal.
      // scale by expected finish time determined by bottleneck
      long flowBitps =
          (*fpIt)->GetBitsLeft() / cf_to_be_schedule->GetAlphaOnline();

      // update utilization profile.
      if (flowBitps > 0) {
        MapWithDef(rates, (*fpIt)->GetFlowId(), flowBitps);
        MapWithInc(sBpsUsed, (*fpIt)->GetSrc(), flowBitps);
        MapWithInc(rBpsUsed, (*fpIt)->GetDest(), flowBitps);
      }
    }
    // Remove capacity from ALL sources and destination for this coflow
    for (map<int, long>::iterator sUsedIt = sBpsUsed.begin();
         sUsedIt != sBpsUsed.end(); sUsedIt++) {
      MapWithDef(sBpsFree, sUsedIt->first, LINK_RATE_BPS);
      sBpsFree[sUsedIt->first] -= sUsedIt->second;
    }
    for (map<int, long>::iterator rUsedIt = rBpsUsed.begin();
         rUsedIt != rBpsUsed.end(); rUsedIt++) {
      MapWithDef(rBpsFree, rUsedIt->first, LINK_RATE_BPS);
      rBpsFree[rUsedIt->first] -= rUsedIt->second;
    }
  }

  // STEP2A: Work conservation as described in the Varys paper.
  RateControlWorkConservationPaper(coflows,
                                   rates,
                                   sBpsFree,
                                   rBpsFree,
                                   LINK_RATE_BPS);
}

// For each ingress port P_in,
// allocate its remaining bandwidth to the coflows in C
// for each active coflow C in Pin, Rem(Pin) is allocated to C’s
//  flows in their current r_ij ratios,
// subject to capacity constraints in corresponding P_out.
void
SchedulerVarysPaper::RateControlWorkConservationPaper(vector<Coflow *> &coflows,
                                                      map<int, long> &rates,
                                                      map<int, long> &sBpsFree,
                                                      map<int, long> &rBpsFree,
                                                      long LINK_RATE_BPS) {

  // map from input port -> total rate allocated to flows behind the input port
  map<int, long> src_rate_map;
  // map from input port -> all flows behind the input port
  map<int, vector<Flow *>> src_flow_map;

  for (vector<Coflow *>::iterator cfIt = coflows.begin();
       cfIt != coflows.end(); cfIt++) {

    if ((*cfIt)->IsComplete()) {
      //such coflow has completed
      continue;
    }

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      long rate_allocated = MapWithDef(rates, (*fpIt)->GetFlowId(), 0L);
      if (rate_allocated > 0) {
        // record rate allocated to this flow
        MapWithInc(src_rate_map, (*fpIt)->GetSrc(), rate_allocated);
        // record this flow after the src
        vector<Flow *> &flows_behind_src = src_flow_map[(*fpIt)->GetSrc()];
        flows_behind_src.push_back(*fpIt);
      }
    }
  }

  for (map<int, vector<Flow *>>::const_iterator kv_iter = src_flow_map.begin();
       kv_iter != src_flow_map.end(); kv_iter++) {
    // For each ingress port P_in,
    int src = kv_iter->first;
    long total_rate_behind_src = MapWithDef(src_rate_map, src, 0L);
    long src_free_bps = MapWithDef(sBpsFree, src, LINK_RATE_BPS);
    const vector<Flow *> &flows_behind_src = kv_iter->second;

    for (vector<Flow *>::const_iterator flow_iter = flows_behind_src.begin();
         flow_iter != flows_behind_src.end(); flow_iter++) {
      // for each active flow after Pin
      // note that flows are in the order of the coflow piority
      // they belong to.
      long previous_rate = MapWithDef(rates, (*flow_iter)->GetFlowId(), 0L);
      long src_rate =
          src_free_bps * (previous_rate / (double) total_rate_behind_src);
      long max_rate_dst =
          MapWithDef(rBpsFree, (*flow_iter)->GetDest(), LINK_RATE_BPS);
      long more_rate = src_rate < max_rate_dst ? src_rate : max_rate_dst;
      if (more_rate > 0) {
        MapWithInc(rates, (*flow_iter)->GetFlowId(), more_rate);
        // update dst free rate on the flight
        MapWithDef(rBpsFree, (*flow_iter)->GetDest(), LINK_RATE_BPS);
        rBpsFree[(*flow_iter)->GetDest()] -= more_rate;
      }
    }
  }

  //debug
  if (DEBUG_LEVEL >= 5) {
    map<int, long> rates_debug = rates;
    cout << " demand " << endl;
    for (vector<Coflow *>::const_iterator cf_iter = coflows.begin();
         cf_iter != coflows.end();
         cf_iter++) {

      Coflow *cf = *cf_iter;
      cout << " coflow id " << cf->GetCoflowId()
           << " alpha " << cf->GetAlpha()
           << " online_alpha " << cf->GetAlphaOnline() << endl;

      vector<Flow *> flows = *(cf->GetFlows());
      for (vector<Flow *>::const_iterator f_iter = flows.begin();
           f_iter != flows.end(); f_iter++) {
        Flow *flow = *f_iter;
        cout << " flow id [" << flow->GetFlowId() << "] "
             << flow->GetSrc() << "->" << flow->GetDest()
             << " demand " << flow->GetBitsLeft()
             << " rate " << MapWithDef(rates_debug, flow->GetFlowId(), (long) 0)
             << endl;
      }
    }
  }
}
