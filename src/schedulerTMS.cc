//
//  SchedulerTMS.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 4/19/15.
//  Copyright (c) 2015 Xin Sunny Huang. All rights reserved.
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


///////////////////////////////////////////////////////
////////////// Code for optc fair based on TMS
///////////////////////////////////////////////////////

SchedulerTMS::SchedulerTMS() : SchedulerOptc() {
}

SchedulerTMS::~SchedulerTMS() {
}

void
SchedulerTMS::Schedule() {

  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerTMS::schedule] scheduling START" << endl;

  Print();

  map<pair<int, int>, long> demandMatrix;

  // STEP 1: Generate demand matrix
  GetDemandAll(demandMatrix);
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);

  // STEP 2: Sinkhorn normalization
  //         Convert the demand matrix into a SQUARE doubly stochastic matrix
  struct timeval sinkhorn_start_time;
  gettimeofday(&sinkhorn_start_time, NULL);
  Sinkhorn(demandMatrix, OPTC_BPS * TMS_CYCLE_LENGTH);
  // Alternative - STEP 2: Stuff the demand into a DSM
  //         Convert the demand matrix into a SQUARE doubly stochastic matrix
  // DemandToDoublySM(demandMatrix);
  struct timeval sinkhorn_end_time;
  gettimeofday(&sinkhorn_end_time, NULL);

  double time_sinkhorn = secondPass(sinkhorn_end_time, sinkhorn_start_time);

  // STEP 3: Generate circuit schedule

  //clearing old patterns
  deepClearSlotVector(m_nextCircuitVector);

  struct timeval decomp_start_time;
  gettimeofday(&decomp_start_time, NULL);
  //generate new patterns
  SchedulerOptc::HopDecomp(demandMatrix, time_sinkhorn);

  struct timeval decomp_end_time;
  gettimeofday(&decomp_end_time, NULL);
  double time_decomp = secondPass(decomp_end_time, decomp_start_time);
  double time_tms = time_sinkhorn + time_decomp;
  if (LOG_COMP_STAT
      && m_coflowPtrVector.size() == 1
      && !m_coflowPtrVector[0]->IsComplete()) {
    // only valid when scheduling for a single coflow.
    m_coflowPtrVector[0]->AddTimeToVector(time_tms);
  }

  if (TMS_LARGEST_N_SLOTS_ONLY > 0) {
    // only apply the largest N slots.
    SortAndSaveLargestNSlots(m_nextCircuitVector, TMS_LARGEST_N_SLOTS_ONLY);
  }

  //debiug
  double total_cycle_length = 0.0;
  for (vector<SLOT *>::const_iterator slot_iter = m_nextCircuitVector.begin();
       slot_iter != m_nextCircuitVector.end(); slot_iter++) {
    total_cycle_length += (*slot_iter)->timeLength;
  }
  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerTMS::schedule] scheduling DONE with "
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
         << "[SchedulerTMS::schedule] new schedule will be activated at "
         << activateTime << "s" << endl;

  } else {

    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerTMS::schedule] no more schedule is available" << endl;
  }
}

//called by schedulerOptc base class
void
SchedulerTMS::CoflowArriveCallBack() {

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
SchedulerTMS::FlowArriveCallBack() {
  // flows should not arrive alone.
}

void
SchedulerTMS::CoflowFinishCallBack(double finishtime) {
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
SchedulerTMS::FlowFinishCallBack(double finishtime) {
  //suspend reschedule upon traffic finish
  //Scheduler::UpdateRescheduleEvent(m_currentTime);

  //however, flow may eat up residual bandwidth
  RateControlOptcFair(m_nextOptcRate, OPTC_BPS);
  // only copy the rate in m_nextElecRate & m_nextOptcRate to flow record.
  // Update flow finish event if needed.
  SetFlowRate();
  // flow finish event will be updated after FlowFinishCallBack().
}

/*
 Sinkhonr in semi-square version
 Perform sinkhorn scaling on DemandMatrix, 
 and reflect the new demand on DemandMatrix in replace 
 such that the col-sum and row-sum equals resulting_sum
 */
void
SchedulerTMS::Sinkhorn(map<pair<int, int>, long> &DemandMatrix,
                       long resulting_sum) {

  map<pair<int, int>, long> origDemand = DemandMatrix;
  map<int, int> src2idx; // src to idx
  map<int, int> dst2idx; // dst to idx
  map<int, int> idx2src; // idx to src
  map<int, int> idx2dst; // idx to dst


  // Step 1: add epslon (small) cells to demand matrix
  //         to make it square
  set<int> srcSet = set<int>();
  set<int> dstSet = set<int>();

  // Step 1.1: find out all src and dst
  for (map<pair<int, int>, long>::iterator dmIt = origDemand.begin();
       dmIt != origDemand.end(); dmIt++) {
    srcSet.insert(dmIt->first.first);
    dstSet.insert(dmIt->first.second);
  }

  // Step 1.2: assign idx to src
  int nodeIdx = 0;
  for (set<int>::iterator srcIt = srcSet.begin();
       srcIt != srcSet.end();
       srcIt++) {
    src2idx[*srcIt] = nodeIdx;
    idx2src[nodeIdx] = *srcIt;
    nodeIdx++;
  }

  // assign idx to dst
  nodeIdx = 0;
  for (set<int>::iterator dstIt = dstSet.begin();
       dstIt != dstSet.end();
       dstIt++) {
    dst2idx[*dstIt] = nodeIdx;
    idx2dst[nodeIdx] = *dstIt;
    nodeIdx++;
  }

  long numSrc = srcSet.size();
  long numDst = dstSet.size();
  long numEnd = max(numSrc, numDst); // num of end point

  double *normDemand = (double *) malloc(sizeof(double) * numEnd * numEnd);

  // Step 1.3: fill with dummy src/dsr for a SQUARE matrix
  // 1/2 either one or none of the branch to be exe to fill demand
  if (numEnd > numSrc) {

    // fill out with (numEnd - numSrc) randomly selected src
    long numSrc2Add = numEnd - numSrc;

    while (numSrc2Add > 0) {
      int randSrc = (rand() % NUM_RACK);
      if (0 > FindWithDef(src2idx, randSrc, -1)) {
        // src available
        int srcIdx = (int) src2idx.size();
        src2idx[randSrc] = srcIdx;
        idx2src[srcIdx] = randSrc;

        //added successfully
        numSrc2Add--;

        //debug
        //cout << "add src " << randSrc << endl;
      }
    }
  }
  // 2/2 either one or none of the branch to be exe to fill demand
  if (numEnd > numDst) {

    // fill out with (numEnd - numSrc) randomly selected dst
    long numDst2Add = numEnd - numDst;

    while (numDst2Add > 0) {
      int randDst = (rand() % NUM_RACK);
      if (0 > FindWithDef(dst2idx, randDst, -1)) {
        // src available
        int dstIdx = (int) dst2idx.size();
        dst2idx[randDst] = dstIdx;
        idx2dst[dstIdx] = randDst;

        //added successfully
        numDst2Add--;

        //debug
        //cout << "add dst " << randDst << endl;
      }
    }
  }

  // convert long demand to double demand
  // and add dummy epslon cells for empty entry
  for (int srcIdx = 0; srcIdx < numEnd; srcIdx++) {
    for (int dstIdx = 0; dstIdx < numEnd; dstIdx++) {
      int src = idx2src[srcIdx];
      int dst = idx2dst[dstIdx];
      *(normDemand + numEnd * srcIdx + dstIdx)
          = (double) FindWithDef(origDemand, make_pair(src, dst), SH_EPSLON);
    }
  }


  // Step 2: sinkhorn normalization
  int iter = 0;
  double err = 1;
  while (iter++ < SH_MAX_ITR_NUM && err > SH_ERROR_BOUND) {

    //norm the rows
    for (int srcIdx = 0; srcIdx < numEnd; srcIdx++) {
      double sum = 0.0;
      for (int dstIdx = 0; dstIdx < numEnd; dstIdx++) {
        sum += *(normDemand + numEnd * srcIdx + dstIdx);
      }
      for (int dstIdx = 0; dstIdx < numEnd; dstIdx++) {
        *(normDemand + numEnd * srcIdx + dstIdx)
            = *(normDemand + numEnd * srcIdx + dstIdx) / sum;
      }
    }
    //norm the cols
    for (int dstIdx = 0; dstIdx < numEnd; dstIdx++) {
      double sum = 0.0;
      for (int srcIdx = 0; srcIdx < numEnd; srcIdx++) {
        sum += *(normDemand + numEnd * srcIdx + dstIdx);
      }
      for (int srcIdx = 0; srcIdx < numEnd; srcIdx++) {
        *(normDemand + numEnd * srcIdx + dstIdx)
            = *(normDemand + numEnd * srcIdx + dstIdx) / sum;
      }
    }

    err = 0.0;
    for (int srcIdx = 0; srcIdx < numEnd; srcIdx++) {
      double sum = 0.0;
      double curnt_error = 0.0;
      for (int dstIdx = 0; dstIdx < numEnd; dstIdx++) {
        sum += *(normDemand + numEnd * srcIdx + dstIdx);
      }
      curnt_error = fabs(1.0 - sum);
      err = max(curnt_error, err);
    }
    for (int dstIdx = 0; dstIdx < numEnd; dstIdx++) {
      double sum = 0.0;
      double curnt_error = 0.0;
      for (int srcIdx = 0; srcIdx < numEnd; srcIdx++) {
        sum += *(normDemand + numEnd * srcIdx + dstIdx);
      }
      curnt_error = fabs(1.0 - sum);
      err = max(curnt_error, err);
    }

  }


  // Step 3: convert double SQUARE doubly-stochastic matrix normDemand
  //         into a long SQUARE doubly-stochastic matrix DemandMatrix
  //         s.t. the demand fits into an interval

  // clear original demand
  DemandMatrix.clear();

  long factor = resulting_sum;
  for (int srcIdx = 0; srcIdx < numEnd; srcIdx++) {
    for (int dstIdx = 0; dstIdx < numEnd; dstIdx++) {
      int src = idx2src[srcIdx];
      int dst = idx2dst[dstIdx];

      double norm = normDemand[srcIdx * numEnd + dstIdx];
      long estDemand = norm * factor;
      MapWithInc(DemandMatrix,
                 make_pair(src, dst),
                 estDemand);

    }
  }
  if (DEBUG_LEVEL >= 1) {
    // sinkhorn ususally has norm error.
    ValidateDSM(DemandMatrix);
  }

  free(normDemand);
}

