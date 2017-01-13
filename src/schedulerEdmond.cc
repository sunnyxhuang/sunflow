//
//  schedulerEdmond.cc
//  coflow-validate-proof
//
//  Created by Xin Sunny Huang on 10/20/15.
//  Copyright © 2015 Xin Sunny Huang. All rights reserved.
//

#include <algorithm>
#include <iomanip>
#include <tgmath.h> // fabs
#include <sys/time.h>

#include "coflow.h"
#include "events.h"
#include "global.h"
#include "scheduler.h"
#include "util.h"


SchedulerEdmond::SchedulerEdmond() : SchedulerOptc() {
}

SchedulerEdmond::~SchedulerEdmond() {
}

void
SchedulerEdmond::Schedule() {

  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerEdmond::schedule] scheduling START" << endl;

  Print();

  struct timeval start_time;
  gettimeofday(&start_time, NULL);

  map<pair<int, int>, long> demandMatrix;

  // STEP 1: Generate demand matrix
  GetDemandAll(demandMatrix);
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);

  struct timeval end_time;
  gettimeofday(&end_time, NULL);

  double demandPrepareTime = secondPass(end_time, start_time);

  // STEP 2: Generate circuit schedule

  //clearing old patterns
  deepClearSlotVector(m_nextCircuitVector);


  //generate a matching based on Edmond's algorithm
  struct timeval edmond_start_time;
  gettimeofday(&edmond_start_time, NULL);

  GenerateOneEdmondAssignment(demandMatrix, demandPrepareTime);

  struct timeval edmond_end_time;
  gettimeofday(&edmond_end_time, NULL);
  double time_edmond = secondPass(edmond_end_time, edmond_start_time);
  if (LOG_COMP_STAT 
      && m_coflowPtrVector.size() == 1
      && !m_coflowPtrVector[0]->IsComplete()) {
    // only valid when scheduling for a single coflow.
    m_coflowPtrVector[0]->AddTimeToVector(time_edmond);
  }

  if (m_nextCircuitVector.size() > 0) {

    // we have new schedule
    m_myTimeLine->RemoveSingularEvent(APPLY_NEW_SCHEDULE);

    // use consumer/producer model
    double activateTime = m_nextCircuitVector[0]->availTime;
    Event *applyScheduleEventPtr = new Event(APPLY_NEW_SCHEDULE, activateTime);
    m_myTimeLine->AddEvent(applyScheduleEventPtr);

    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerEdmond::schedule] new schedule will be activated at "
         << activateTime << "s" << endl;

  } else {

    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerEdmond::schedule] no more schedule is available" << endl;
  }
}

void
SchedulerEdmond::GenerateOneEdmondAssignment(map<pair<int, int>, long> &demandMatrix,
                                             double timeSinceCalStart) {
  if (!m_nextCircuitVector.empty() && !demandMatrix.empty()) {
    cout << "try to generate new assignments "
         << " for non-empty next circuit vector!" << endl;
  }

  cout << "now generating patterns" << endl;

  struct timeval pattern_cal_start_time;
  struct timeval pattern_cal_end_time;

  set<pair<int, int>> *circuits = new set<pair<int, int>>();

  // STEP 1: get start time
  gettimeofday(&pattern_cal_start_time, NULL);

  // STEP 2: find a matching
  //******* Edmond matching begin *********
  //                                          *
  DemandMatchingEdmond(demandMatrix, circuits);
  //                                          *
  //******* Edmond matching end ***********

  gettimeofday(&pattern_cal_end_time, NULL);

  // STEP 3: find circuit assignment metrics

  // STEP 3.2: slot length based on demand served
  double slot_length = EDMOND_CYCLE_LENGTH;

  // STEP 3.3: slot efficiency score
  double slot_score = GetSlotEfficiency(circuits, demandMatrix);

  // STEP 3.4: control plan delay
  double patternCalTime = secondPass(pattern_cal_end_time, pattern_cal_start_time);

  // test of pattern computation delay
  if (ZERO_COMP_TIME) {
    timeSinceCalStart = 0.0;
    patternCalTime = 0.0;
  }

  timeSinceCalStart += patternCalTime;
  double pattern_avail_time = m_currentTime + timeSinceCalStart;

  // STEP 3.5: a new assignment is generated!
  SLOT *slotp = new SLOT(slot_length, circuits, slot_score, pattern_avail_time);
  m_nextCircuitVector.push_back(slotp);

}

//called by schedulerOptc base class
void
SchedulerEdmond::CoflowArriveCallBack() {

  if (m_nextCircuitVector.empty() && m_circuitVector.empty()) {
    // kick off scheduling.
    Scheduler::UpdateRescheduleEvent(m_currentTime);
    return;
  }

  // otherwise, we are in the middle of scheduling
  //  suspend reschedule upon traffic arrival

  //however, flow may take up a share of bandwidth
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);
  // only copy the rate in m_nextElecRate & m_nextOptcRate to flow record.
  // Update flow finish event if needed.
  SetFlowRate();
  UpdateFlowFinishEvent(m_currentTime);
}

//wrapper called by traffic generator

void
SchedulerEdmond::FlowArriveCallBack() {
}

void
SchedulerEdmond::CoflowFinishCallBack(double finishtime) {
  //suspend reschedule upon traffic finish
  //Scheduler::UpdateRescheduleEvent(m_currentTime);

  //however, flow may eat up residual bandwidth
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);
  // only copy the rate in m_nextElecRate & m_nextOptcRate to flow record.
  // Update flow finish event if needed.
  SetFlowRate();
  // flow finish event will be updated after FlowFinishCallBack().
}

void
SchedulerEdmond::FlowFinishCallBack(double finishtime) {
  //suspend reschedule upon traffic finish
  //Scheduler::UpdateRescheduleEvent(m_currentTime);

  //however, flow may eat up residual bandwidth
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);
  // only copy the rate in m_nextElecRate & m_nextOptcRate to flow record.
  // Update flow finish event if needed.
  SetFlowRate();
  // flow finish event will be updated after FlowFinishCallBack().
}

