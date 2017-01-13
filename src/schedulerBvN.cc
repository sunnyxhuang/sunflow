//
//  schedulerBvN.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 12/14/15.
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
////////////// Code for optc fair based on BvN (optimal under Electrical)
///////////////////////////////////////////////////////////////////

SchedulerBvN::SchedulerBvN() : SchedulerOptc() {
}

SchedulerBvN::~SchedulerBvN() {
}

void
SchedulerBvN::Schedule() {

  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerBvN::schedule] scheduling START" << endl;

  Print();

  struct timeval start_time;
  gettimeofday(&start_time, NULL);

  map<pair<int, int>, long> demandMatrix;

  // STEP 1: Generate demand matrix
  GetDemandAll(demandMatrix);
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);

  // STEP 2: Stuff the demand into a DSM
  //         Convert the demand matrix into a SQUARE doubly stochastic matrix
  DemandToDoublySM(demandMatrix);

  struct timeval end_time;
  gettimeofday(&end_time, NULL);

  double demandPrepareTime = secondPass(end_time, start_time);

  // STEP 3: Generate circuit schedule

  //clearing old patterns
  deepClearSlotVector(m_nextCircuitVector);

  //generate new patterns
  SchedulerOptc::HopDecomp(demandMatrix, demandPrepareTime);

  //debiug
  double total_cycle_length = 0.0;
  for (vector<SLOT *>::const_iterator slot_iter = m_nextCircuitVector.begin();
       slot_iter != m_nextCircuitVector.end(); slot_iter++) {
    total_cycle_length += (*slot_iter)->timeLength;
  }
  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerBvN::schedule] scheduling DONE with "
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
         << "[SchedulerBvN::schedule] new schedule will be activated at "
         << activateTime << "s" << endl;

  } else {

    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerBvN::schedule] no more schedule is available" << endl;
  }
}

//called by schedulerOptc base class
void
SchedulerBvN::CoflowArriveCallBack() {
  Scheduler::UpdateRescheduleEvent(m_currentTime);
}

//wrapper called by traffic generator

void
SchedulerBvN::FlowArriveCallBack() {
  // flows should not arrive alone.
}

void
SchedulerBvN::CoflowFinishCallBack(double finishtime) {
  Scheduler::UpdateRescheduleEvent(finishtime);
}

void
SchedulerBvN::FlowFinishCallBack(double finishtime) {
  //suspend reschedule upon traffic finish
  //Scheduler::UpdateRescheduleEvent(m_currentTime);

  //however, flow may eat up residual bandwidth
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);
  // only copy the rate in m_nextElecRate & m_nextOptcRate to flow record.
  // Update flow finish event if needed.
  SetFlowRate();
  // flow finish event will be updated after FlowFinishCallBack().
}
