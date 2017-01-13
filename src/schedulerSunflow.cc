//
//  schedulerSunflow.cpp
//  Ximulator
//
//  Created by Xin Sunny Huang on 11/5/15.
//  Copyright © 2015 Xin Sunny Huang. All rights reserved.
//

#include <algorithm>
#include <iomanip>
#include <sys/time.h>

#include "coflow.h"
#include "events.h"
#include "global.h"
#include "scheduler.h"
#include "sunflowSolver.h"
#include "util.h"

///////////////////////////////////////////////////////
////////////// Code for Sunflow
///////////////////////////////////////////////////////


SchedulerSunflow::SchedulerSunflow() : SchedulerOptc() {
}

SchedulerSunflow::~SchedulerSunflow() {
}

void
SchedulerSunflow::ApplyNewSchedule() {
  if (m_reservations.empty()) {
    cout
        << "[SchedulerSunflow::ApplyNewSchedule] error: No new schedule is available"
        << endl;
    return;
  }

  circuits_to_apply_list.clear();
  m_circuits_switching_to.clear();

  // reset expiration time upon new schedule.
  m_circuit_to_expired_time.clear();

  // set flow rate based on m_nextElecRate & m_nextOptcRate.
  // flow finish event will be updated in OrderCircuit();
  SetFlowRate();

  // only one ORDER_CIRCUIT event at a time.
  m_myTimeLine->RemoveSingularEvent(ORDER_CIRCUIT);

  // there could be multiple APPLY_CIRCUIT events.
  m_myTimeLine->RemoveMultipleEvent(APPLY_CIRCUIT);

  // update schedule start time
  // this_schedule_start_time = m_currentTime;

  // update schedule end event
  m_myTimeLine->RemoveSingularEvent(SCHEDULE_END);
  double scheduleEndTime = m_currentTime + m_last_schedule_valid_for;
  Event *scheduleEndEventPtr = new Event(SCHEDULE_END, scheduleEndTime);
  m_myTimeLine->AddEvent(scheduleEndEventPtr);

  OrderCircuit();
}

bool
SchedulerSunflow::OrderCircuit(void) {
  if (m_reservations.empty()) {
    cout << "[SchedulerSunflow::OrderCircuit]"
         << " m_port_reserve is empty => no scheduling" << endl;
    return false;
  }

  //update circuits
  circuits_to_apply_list.push_back(std::make_pair(
      m_currentTime + CIRCUIT_RECONF_DELAY,
      set<pair<int, int>>()));
  set<pair<int, int>> &circuits_to_set = circuits_to_apply_list.back().second;

  double this_reservation_virtual_time =
      GetNextReservationTimeAndCircuits(circuits_to_set);

  set<pair<int, int>> remaining_circuits;
  ExcludeConflictingCircuits(m_activeCircuits,
                             circuits_to_set,
                             remaining_circuits);

  m_activeCircuits.clear();
  m_activeCircuits = remaining_circuits;

  SetFlowPath();
  UpdateFlowFinishEvent(m_currentTime);

  // update switching circuits set.
  m_circuits_switching_to.insert(circuits_to_set.begin(),
                                 circuits_to_set.end());

  // add apply circuit event.
  // APPLY_CIRCUIT event should be added BEFORE next reservation alarm is added.
  Event *applyCircuitEventPtr = new Event(APPLY_CIRCUIT,
                                          m_currentTime + CIRCUIT_RECONF_DELAY);
  m_myTimeLine->AddEvent(applyCircuitEventPtr);

  // add next reservation event, if any.
  if (!m_reservations.empty()) {
    double next_reservation_virtual_time = m_reservations.top().first;
    double next_reservation_time = m_currentTime
        + (next_reservation_virtual_time - this_reservation_virtual_time);
    Event
        *orderCircuitEventPtr = new Event(ORDER_CIRCUIT, next_reservation_time);
    m_myTimeLine->AddEvent(orderCircuitEventPtr);
  }

  if (DEBUG_LEVEL >= 9) {
    cout << __func__ << " active circuits "
         << CircuitsToString(m_activeCircuits) << endl;
  }

  // debug
  // check for resource leak.
  if (ENABLE_LEAK_CHECK) {
    set<pair<int, int>> next_active_circuits = m_activeCircuits;
    next_active_circuits.insert(circuits_to_set.begin(), circuits_to_set.end());
    double tolerance = GetLeakCheckTolerance();
    if (!ValidateNoResourceLeak(m_coflowPtrVector,
                                next_active_circuits,
                                tolerance)) {
      cout << "Error in " << __func__ << ": See resource leakage!!!" << endl;
      if (LEAK_CHECK_EXIT) {
        exit(-1);
      }
    }
  }

  return true; // circuits changed
}

double
SchedulerSunflow::GetNextReservationTimeAndCircuits(set<pair<int,
                                                             int>> &circuits_to_set) {
  circuits_to_set.clear();
  if (m_reservations.empty()) {
    return 0;
  }
  double next_time = m_reservations.top().first;
  while (!m_reservations.empty() && next_time >= m_reservations.top().first) {
    next_time = m_reservations.top().first;
    for (map<pair<int, int>, double>::const_iterator
             circuit_time_pair = m_reservations.top().second.begin();
         circuit_time_pair != m_reservations.top().second.end();
         circuit_time_pair++) {
      circuits_to_set.insert(circuit_time_pair->first);

      // update expiration time when new circuits
      // are obtained from m_reservations
      // and ready to deploy.
      // expire_time should base on current time
      // and time in m_reservations is virtual time.
      // NOTE : MAKE SURE this function is called when
      // real m_currentTime == virtual next_time
      double
          expire_time = circuit_time_pair->second - next_time + m_currentTime;
      m_circuit_to_expired_time[circuit_time_pair->first] = expire_time;
    }
    m_reservations.pop();
  }
  return next_time;
}

// called upon flow finish event to mark potential leakage.
void
SchedulerSunflow::MarkPotentialResourceLeak(double alarm_time) {
  set<pair<int, int>> potential_leaked_circuits;
  FindLeakedCircuits(m_coflowPtrVector,
                     m_activeCircuits,
                     potential_leaked_circuits);
  if (!potential_leaked_circuits.empty()) {
    m_leaked_circuits_profile.push_back(std::make_pair(alarm_time,
                                                       potential_leaked_circuits));
  }
}

// find circuits with demand
// whose src/dst are not taken up.
void
SchedulerSunflow::FindLeakedCircuits(const vector<Coflow *> &coflows,
                                     const set<pair<int,
                                                    int>> &occupied_circuits,
                                     set<pair<int, int>> &leaked_circuits) {

  set<pair<int, int>> circuits_with_demand;
  FindCircuitsWithDemand(coflows, circuits_with_demand);

  set<int> busy_src;
  set<int> busy_dst;
  for (set<pair<int, int>>::const_iterator
           occupied_circuit = occupied_circuits.begin();
       occupied_circuit != occupied_circuits.end();
       occupied_circuit++) {
    if (ContainsKey(circuits_with_demand, *occupied_circuit)) {
      busy_src.insert(occupied_circuit->first);
      busy_dst.insert(occupied_circuit->second);
    }
  }

  for (set<pair<int, int>>::const_iterator
           demanding_circuit = circuits_with_demand.begin();
       demanding_circuit != circuits_with_demand.end();
       demanding_circuit++) {
    if (!ContainsKey(busy_src, demanding_circuit->first)
        && !ContainsKey(busy_dst, demanding_circuit->second)) {
      leaked_circuits.insert(*demanding_circuit);
    }
  }

}

void
SchedulerSunflow::FindCircuitsWithDemand(const vector<Coflow *> &coflows,
                                         set<pair<int,
                                                  int>> &circuits_with_demand) {
  for (vector<Coflow *>::const_iterator cfIt = coflows.begin();
       cfIt != coflows.end(); cfIt++) {
    if ((*cfIt)->IsComplete() || 0 == (*cfIt)->GetAlpha()) {
      // such coflow has completed
      // or has zero demand
      continue;
    }

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();
    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }
      int src = (*fpIt)->GetSrc();
      int dst = (*fpIt)->GetDest();

      circuits_with_demand.insert(std::make_pair(src, dst));
    }
  }
}

double
SchedulerSunflow::GetLeakCheckTolerance() {
  double tolerance = CIRCUIT_RECONF_DELAY + 1000 / (double) OPTC_BPS;
  if (DEBUG_LEVEL >= 30) {
    // use a smaller tolerance for debugging.
    tolerance = CIRCUIT_RECONF_DELAY;
  }
  if (tolerance < 0.01) {
    // never use 0 tolerance to avoid false alarm.
    tolerance = 0.01;
  }
  return tolerance;
}

// called in order_circuit event to validate leakage.
bool
SchedulerSunflow::ValidateNoResourceLeak(
    const vector<Coflow *> &coflows,
    const set<pair<int, int>> &next_active_circuits,
    double tolerance) {
  // TODO: implement resource leakage check.
  // It now has some logical problem.
  // return true;
  bool validate_ok = true;
  // check for leakage.
  for (vector<pair<double, set<pair<int, int>>>>::iterator
           potential_leaked_profile = m_leaked_circuits_profile.begin();
       potential_leaked_profile != m_leaked_circuits_profile.end();
       potential_leaked_profile++) {
    double leaked_since = potential_leaked_profile->first;
    if (leaked_since + tolerance < m_currentTime) {
      cout << "circuits have been idle during (" << leaked_since << ", "
           << m_currentTime << ") : ";
      set<pair<int, int>> &leaked_circuits = potential_leaked_profile->second;
      for (set<pair<int, int>>::const_iterator
               circuit = leaked_circuits.begin();
           circuit != leaked_circuits.end();
           circuit++) {
        cout << circuit->first << "->" << circuit->second << ", ";
      }
      cout << endl;
      validate_ok = false;
    }
  }

  // erase history if extra circuit are no longer valid.
  // not valid => no demand or conflicting with current active circuits
  for (vector<pair<double, set<pair<int, int>>>>::iterator
           potential_leaked_profile = m_leaked_circuits_profile.begin();
       potential_leaked_profile != m_leaked_circuits_profile.end();
    /*step advance missing*/) {

    set<pair<int, int>>
        &leaked_circuit_history = potential_leaked_profile->second;
    // not_conflicting_leaked_circuits = leaked_circuit_history - current_active_circuits
    set<pair<int, int>> not_conflicting_leaked_circuits;
    // cout << __func__ << " calling ExcludeConflictingCircuits " << endl;
    ExcludeConflictingCircuits(leaked_circuit_history,
                               next_active_circuits,
                               not_conflicting_leaked_circuits);
    // remaining_leaked_circuits = not_conflicting_leaked_circuits - circuits_with_demand
    set<pair<int, int>> circuits_with_demand;
    FindCircuitsWithDemand(coflows, circuits_with_demand);
    // cout << __func__ << " circuits with demand "
    //    << CircuitsToString(circuits_with_demand) << endl;
    set<pair<int, int>> remaining_leaked_circuits;
    SetMinus(not_conflicting_leaked_circuits,
             circuits_with_demand,
             remaining_leaked_circuits);

    if (remaining_leaked_circuits.empty()) {
      potential_leaked_profile =
          m_leaked_circuits_profile.erase(potential_leaked_profile);
    } else {
      potential_leaked_profile->second.swap(remaining_leaked_circuits);
      potential_leaked_profile++;
    }
  }

  return validate_ok;
}

// remove circuits from current_circuits
// that that not in circuits_to_set,
// but are using src/dst requested in circuits_to_set.
void
SchedulerSunflow::ExcludeConflictingCircuits(
    const set<pair<int, int> > &current_circuits,
    const set<pair<int, int> > &circuits_to_set,
    set<pair<int, int> > &remaining_circuits) {

  set<int> excluding_src;
  set<int> excluding_dst;
  for (set<pair<int, int>>::const_iterator
           new_circuit = circuits_to_set.begin();
       new_circuit != circuits_to_set.end();
       new_circuit++) {

    excluding_src.insert(new_circuit->first);
    excluding_dst.insert(new_circuit->second);
  }

  remaining_circuits.clear();
  // initialize remaining_circuits to current_circuits.
  remaining_circuits = current_circuits;

  for (set<pair<int, int>>::const_iterator
           current_circuit = current_circuits.begin();
       current_circuit != current_circuits.end();
       current_circuit++) {
    if (!ContainsKey(circuits_to_set, *current_circuit)
        // current circuit not in new circuits.
        && (ContainsKey(excluding_src, current_circuit->first)
            || ContainsKey(excluding_dst, current_circuit->second))) {
      // and current circuit is in conflict.)
      remaining_circuits.erase(*current_circuit);
    }
  }

  if (DEBUG_LEVEL >= 15) {
    cout << __func__ << endl;
    cout << "   current circuits " << CircuitsToString(current_circuits)
         << endl;
    cout << "   new circuits " << CircuitsToString(circuits_to_set) << endl;
    cout << "   remaining circuits " << CircuitsToString(remaining_circuits)
         << endl;
  }

}

// lhs_minus_rhs = circuits in lhs, but not in rhs.
void
SchedulerSunflow::SetMinus(const set<pair<int, int>> &lhs,
                           const set<pair<int, int>> &rhs,
                           set<pair<int, int>> &lhs_minus_rhs) {
  for (set<pair<int, int>>::const_iterator
           circuit_in_lhs = lhs.begin();
       circuit_in_lhs != lhs.end();
       circuit_in_lhs++) {
    if (!ContainsKey(rhs, *circuit_in_lhs)) {
      lhs_minus_rhs.insert(*circuit_in_lhs);
    }
  }
}

bool
SchedulerSunflow::ApplyCircuit(void) {
  if (circuits_to_apply_list.empty()) {
    cout << "circuits_to_apply_list is empty while applying new circuits"
         << endl;
    return false;
  }

  if (m_currentTime < circuits_to_apply_list.front().first) {
    cout << "error: inconsistent circuit apply time " << endl;
    return false;
  }

  // m_currentTime > circuit valid time

  set<pair<int, int>>
      &circuits_newly_added = circuits_to_apply_list.front().second;
  // debug
  if (DEBUG_LEVEL >= 15) {
    cout << m_currentTime << "s adding circuits "
         << CircuitsToString(circuits_newly_added) << endl;
  }

  set<pair<int, int>> next_active_circuits = m_activeCircuits;
  next_active_circuits.insert(circuits_newly_added.begin(),
                              circuits_newly_added.end());

  bool is_valid_circuits = ValidateCircuitAssignment(next_active_circuits);
  if (!is_valid_circuits) {
    cout << "Error: Circuit assignment is not valid" << endl;
    exit(-1);
  }

  // debug
  // check for resource leak.
  if (ENABLE_LEAK_CHECK) {
    double tolerance = GetLeakCheckTolerance();
    if (!ValidateNoResourceLeak(m_coflowPtrVector,
                                next_active_circuits,
                                tolerance)) {
      cout << "Error in " << __func__ << ": See resource leakage!!!" << endl;
      if (LEAK_CHECK_EXIT) {
        exit(-1);
      }
    }
  }

  // update active circuits.
  m_activeCircuits.swap(next_active_circuits);

  circuits_to_apply_list.erase(circuits_to_apply_list.begin());

  SetFlowPath();
  UpdateFlowFinishEvent(m_currentTime);

  // update switching circuits set.
  // still_switching_circuits = m_circuits_switching_to - m_activeCircuits
  set<pair<int, int>> still_switching_circuits;
  for (set<pair<int, int>>::const_iterator
           switching_circuit = m_circuits_switching_to.begin();
       switching_circuit != m_circuits_switching_to.end();
       switching_circuit++) {
    if (!ContainsKey(m_activeCircuits, *switching_circuit)) {
      // this switching circuit is still switching (NOT active yet).
      still_switching_circuits.insert(*switching_circuit);
    }
  }
  m_circuits_switching_to.clear();
  m_circuits_switching_to = still_switching_circuits;

  if (DEBUG_LEVEL >= 9) {
    cout << __func__ << " active circuits "
         << CircuitsToString(m_activeCircuits) << endl;
  }

  return true; // circuits changed
}

void
SchedulerSunflow::Schedule() {
  // debug
  cout << fixed << setw(FLOAT_TIME_WIDTH) << m_currentTime << "s "
       << "[SchedulerSunflow::schedule] sunflow starts. " << endl;
  cout << " Coflow profiles " << endl;
  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {
    cout << "            " << (*cfIt)->toString() << endl;
  }

  // Step 1: each coflow reserves circuits.
  struct timeval start_time;
  gettimeofday(&start_time, NULL);

  ComputeReservations(m_coflowPtrVector);

  struct timeval end_time;
  gettimeofday(&end_time, NULL);
  double prepare_time = secondPass(end_time, start_time);

  if (ZERO_COMP_TIME) prepare_time = 0.0;

  // Step 3: perform rate control.
  SchedulerOptc::RateControlBySortedCoflows(m_coflowPtrVector, m_nextOptcRate);

  // Step 4: update new pattern event if necessary.
  if (!m_reservations.empty()) {

    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerSunflow::schedule] sunflow done with "
         << setw(5) << m_reservations.size() << " circuit sets / "
         << m_last_schedule_valid_for << " seconds" << endl;

    if (SPEEDUP_1BY1_TRAFFIC_IN_SCHEDULER
        && m_coflowPtrVector.size() == 1) {
      // sunny is hacking!
      vector<Flow *> no_finished_flows;
      ScheduleToNotifyTrafficFinish(m_currentTime + m_last_schedule_valid_for,
                                    m_coflowPtrVector,
                                    no_finished_flows);
      m_coflowPtrVector.clear();
      return;
    }

    // we have new schedule
    m_myTimeLine->RemoveSingularEvent(APPLY_NEW_SCHEDULE);
    double activateTime = m_currentTime + prepare_time;
    Event *applyScheduleEventPtr = new Event(APPLY_NEW_SCHEDULE, activateTime);
    m_myTimeLine->AddEvent(applyScheduleEventPtr);

    // maintain current active circuits for better connectivity.

  } else {

    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerSunflow::schedule] no more schedule is available"
         << endl;
  }
}

void
SchedulerSunflow::ComputeReservations(vector<Coflow *> &coflows) {
  m_last_schedule_valid_for = 0;
  // clear all previous reservations.
  while (!m_reservations.empty()) {
    m_reservations.pop();
  }

  if (coflows.empty()) return;

  // convert coflow demand to layered
  vector<map<pair<int, int>, long>> layered_demand;
  CoflowsToLayeredDemands(coflows, layered_demand);
  if (layered_demand.empty()) return;

  SunflowSolver sunflow_solver;
  sunflow_solver.ReserveCircuits(layered_demand);
  m_last_schedule_valid_for = sunflow_solver.GetMaxReservationTime();

  if (LOG_COMP_STAT
      && coflows.size() == 1) {
    CompTimeBreakdown comptime;
    comptime.m_time_reserve = sunflow_solver.GetCompTimeReserve();
    coflows[0]->SetCompTime(comptime);
  }

  map<pair<int, int>, double> circuits_to_expire_time;
  double next_reservation_time;
  while (sunflow_solver.HasNext()) {
    circuits_to_expire_time.clear();
    next_reservation_time
        = sunflow_solver.GetNextCircuitsAndExpireTime(circuits_to_expire_time);
    m_reservations.push(std::make_pair(next_reservation_time,
                                       circuits_to_expire_time));
  }
}

//called by schedulerOptc base class
void
SchedulerSunflow::CoflowArriveCallBack() {
  Scheduler::UpdateRescheduleEvent(m_currentTime);
  // sorts coflows upon arrival and departure.
  CalAlphaAndSortCoflowsInPlace(m_coflowPtrVector);
}

void
SchedulerSunflow::FlowArriveCallBack() {
  // flows should not arrive alone.
}

void
SchedulerSunflow::CoflowFinishCallBack(double finishtime) {
  Scheduler::UpdateRescheduleEvent(finishtime);
  // only sorts coflows upon arrival and departure.
  CalAlphaAndSortCoflowsInPlace(m_coflowPtrVector);
}

void
SchedulerSunflow::FlowFinishCallBack(double finishtime) {
  // Rate controlled is performed in the rescheduling
  SchedulerOptc::RateControlBySortedCoflows(m_coflowPtrVector,
                                            m_nextOptcRate);
  SetFlowRate();
  // flow finish event will be updated in SchedulerOptc::Transmit()
  // after flow is finished.
  // UpdateFlowFinishEvent(finishtime);

  //debug
  if (ENABLE_LEAK_CHECK) {
    MarkPotentialResourceLeak(finishtime);
  }
}
