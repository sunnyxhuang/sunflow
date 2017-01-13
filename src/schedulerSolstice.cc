//
//  schedulerSolstice.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 12/11/15.
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


///////////////////////////////////////////////////////////////////
////////////// Code for Solstice (coNext'15)
///////////////////////////////////////////////////////////////////

SchedulerSolstice::SchedulerSolstice() : SchedulerOptc() {
}

SchedulerSolstice::~SchedulerSolstice() {
}

void
SchedulerSolstice::Schedule() {

  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerSolstice::schedule] scheduling START" << endl;

  Print();

  // STEP 1: Generate demand matrix
  map<pair<int, int>, long> demandMatrix;
  GetDemandAll(demandMatrix);
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);

  // STEP 2: Stuff the demand into a DSM
  //         Convert the demand matrix into a SQUARE doubly stochastic matrixs
  struct timeval stuff_start_time;
  gettimeofday(&stuff_start_time, NULL);
  // "quick stuff"
  DemandToDoublySM(demandMatrix);

  struct timeval stuff_end_time;
  gettimeofday(&stuff_end_time, NULL);
  double time_stuff = secondPass(stuff_end_time, stuff_start_time);

  // STEP 3: Generate circuit schedule

  //clearing old patterns
  deepClearSlotVector(m_nextCircuitVector);

  //generate new patterns
  struct timeval slice_start_time;
  gettimeofday(&slice_start_time, NULL);
  // do valid work
  SchedulerOptc::BigSliceDecomp(demandMatrix, time_stuff);
  struct timeval slice_end_time;
  gettimeofday(&slice_end_time, NULL);
  double time_slice = secondPass(slice_end_time, slice_start_time);
  if (LOG_COMP_STAT
      && m_coflowPtrVector.size() == 1
      && !m_coflowPtrVector[0]->IsComplete()) {
    // only valid when scheduling for a single coflow.
    CompTimeBreakdown compstat;
    compstat.m_time_stuff = time_stuff;
    compstat.m_time_slice = time_slice;
    m_coflowPtrVector[0]->SetCompTime(compstat);
  }

  //debug
  double total_cycle_length = 0.0;
  for (vector<SLOT *>::const_iterator slot_iter = m_nextCircuitVector.begin();
       slot_iter != m_nextCircuitVector.end(); slot_iter++) {
    total_cycle_length += (*slot_iter)->timeLength;
  }
  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerSolstice::schedule] scheduling DONE with "
       << m_nextCircuitVector.size() << " slots, "
       << total_cycle_length << "s in total" << endl;

  if (m_nextCircuitVector.size() > 0) {

    // we have new schedule
    m_myTimeLine->RemoveSingularEvent(APPLY_NEW_SCHEDULE);

    // use consumer/producer model
    double activateTime = m_nextCircuitVector[0]->availTime;
    Event *applyScheduleEventPtr = new Event(APPLY_NEW_SCHEDULE, activateTime);
    m_myTimeLine->AddEvent(applyScheduleEventPtr);

    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerSolstice::schedule] new schedule will be activated at "
         << activateTime << "s" << endl;

  } else {

    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerSolstice::schedule] no more schedule is available" << endl;
  }
}

//called by schedulerOptc base class
void
SchedulerSolstice::CoflowArriveCallBack() {
  Scheduler::UpdateRescheduleEvent(m_currentTime);
}

//wrapper called by traffic generator

void
SchedulerSolstice::FlowArriveCallBack() {
  // flows should not arrive alone.
}

void
SchedulerSolstice::CoflowFinishCallBack(double finishtime) {
  Scheduler::UpdateRescheduleEvent(finishtime);
}

void
SchedulerSolstice::FlowFinishCallBack(double finishtime) {
  //suspend reschedule upon traffic finish
  //Scheduler::UpdateRescheduleEvent(m_currentTime);

  //however, flow may eat up residual bandwidth
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);
  // only copy the rate in m_nextElecRate & m_nextOptcRate to flow record.
  // Update flow finish event if needed.
  SetFlowRate();
  // flow finish event will be updated after FlowFinishCallBack().
}


