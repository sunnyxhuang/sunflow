//
//  schedulerAaloImpl.cc
//  Ximulator
//
//  Reimplementation of Aalo as seen on Github.
//

#include <algorithm>
#include <iomanip>
#include <sys/time.h>

#include "scheduler.h"
#include "events.h"
#include "global.h"
#include "util.h"
#include "coflow.h"

SchedulerAaloImpl::SchedulerAaloImpl() : SchedulerVarys() {
  m_coflow_jid_queues.resize(AALO_Q_NUM);
}

void
SchedulerAaloImpl::AddCoflows(vector<Coflow *> *cfsPtr) {
  if (!cfsPtr) return;
  // check all coflows are legal
  for (vector<Coflow *>::iterator cfpIt = cfsPtr->begin();
       cfpIt != cfsPtr->end(); cfpIt++) {
    if (!(*cfpIt)->isRawCoflow()) {
      cout << "[SchedulerAaloImpl::AddCoflows]"
           << "error: is not a raw coflow" << endl;
      continue;
    }
    m_coflowPtrVector.push_back(*cfpIt);
    // add to the highest piority queue.
    m_coflow_jid_queues[0].push_back((*cfpIt)->GetJobId());
  }
}

void
SchedulerAaloImpl::Schedule() {
  cout << fixed << setw(FLOAT_TIME_WIDTH)
       << m_currentTime << "s "
       << "[SchedulerAaloImpl::schedule] AaloImpl scheduling START" << endl;

  Print();

  struct timeval start_time;
  gettimeofday(&start_time, NULL);
  /////////////// Aalo //////////////////////////

  // STEP 1: Initialize next rate for all flows to (0,0)
  m_nextElecRate.clear();
  m_nextOptcRate.clear();

  // STEP 2: Perform Aalo rate control
  RateControlAaloImpl(m_coflowPtrVector, m_coflow_jid_queues,
                      m_nextElecRate, ELEC_BPS);

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

// perform aalo based rate control - as seen in Github.
// "FIFO within each queue Strict priority across queues" - Github Aalo impl
// update Coflow queue/piority based on data sent on the Coflow
// perform max-min among flows within the same queue
// recored flow rate to rates (in place).
// parameters@
// coflows: a list of coflows
// rates : map from flow id to allocated rate.
void
SchedulerAaloImpl::RateControlAaloImpl(vector<Coflow *> &coflows,
                                       vector<vector<int>> &coflow_id_queues,
                                       map<int, long> &rates,
                                       long LINK_RATE_BPS) {
  if (coflows.empty()) {
    return;
  }

  rates.clear();

  map<int, Coflow *> coflow_id_ptr_map;
  // update Coflow piority/queue.
  // coflow are sorted by
  UpdateCoflowQueue(coflows, coflow_id_queues, coflow_id_ptr_map);

  // initialize.
  map<int, long> sBpsFree, rBpsFree;
  // for each queue
  for (unsigned int queue_runner = 0;
       queue_runner < AALO_Q_NUM;
       queue_runner++) {
    vector<int> &current_queue = coflow_id_queues[queue_runner];

    map<int, int> src_flow_num;
    map<int, int> dst_flow_num;

    // for each coflow
    for (vector<int>::iterator
             coflow_in_q = current_queue.begin();
         coflow_in_q != current_queue.end();
         coflow_in_q++) {

      Coflow *cfp = coflow_id_ptr_map[*coflow_in_q];

      if (cfp->IsFlowsAddedComplete()) {
        // coflow not completed yet
        // but the flows added so far have finished
        continue;
      }

      map<int, long> sBpsUsed, rBpsUsed;

      vector<Flow *> *flowVecPtr = cfp->GetFlows();
      // calculate flow number on sender and receiver side.
      // for each flow within coflow
      for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
           fpIt != flowVecPtr->end(); fpIt++) {

        if ((*fpIt)->GetBitsLeft() <= 0) {
          // such flow has finished
          continue;
        }

        int src = (*fpIt)->GetSrc();
        int dst = (*fpIt)->GetDest();
        MapWithInc(src_flow_num, src, 1);
        MapWithInc(dst_flow_num, dst, 1);
      } // for each flow within coflow

      // Determine rate based only on this job and available bandwidth
      // for each flow within coflow
      for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
           fpIt != flowVecPtr->end(); fpIt++) {

        if ((*fpIt)->GetBitsLeft() <= 0) {
          // such flow has finished
          continue; // with next flow
        }

        int src = (*fpIt)->GetSrc();
        int dst = (*fpIt)->GetDest();

        long sBps = MapWithDef(sBpsFree, src, LINK_RATE_BPS);
        long rBps = MapWithDef(rBpsFree, dst, LINK_RATE_BPS);

        if (sBps <= 0 || rBps <= 0) {
          // no more bandwidth left.
          continue; // with next flow
        }

        long srcAvgBps = sBps / src_flow_num[src];
        long dstAvgBps = rBps / dst_flow_num[dst];
        long flowBitps = min(srcAvgBps, dstAvgBps);

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
    } // for each coflow

  } // for each queue


  //debug
  if (DEBUG_LEVEL >= 5) {
    map<int, long> rates_debug = rates;
    cout << " demand " << endl;
    for (vector<Coflow *>::const_iterator
             cf_iter = coflows.begin();
         cf_iter != coflows.end();
         cf_iter++) {

      Coflow *cf = *cf_iter;
      cout << " coflow Jid " << cf->GetJobId()
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
// build up coflow_id_ptr_map for later reference
// clear obsolete coflow from last_coflow_id_queues
//
void
SchedulerAaloImpl::UpdateCoflowQueue(
    vector<Coflow *> &coflows,
    vector<vector<int>> &coflow_id_queues,
    map<int, Coflow *> &active_coflow_id_ptr_map) {

  active_coflow_id_ptr_map.clear();

  // coflow_id_ptr_map has all active coflows
  for (vector<Coflow *>::const_iterator
           coflow_iter = coflows.begin();
       coflow_iter != coflows.end();
       coflow_iter++) {
    int coflow_jid = (*coflow_iter)->GetJobId();

    active_coflow_id_ptr_map[coflow_jid] = (*coflow_iter);
  }

  // update the queue
  for (unsigned int queue_runner = 0;
       queue_runner < AALO_Q_NUM;
       queue_runner++) {
    vector<int> &current_queue = coflow_id_queues[queue_runner];
    for (vector<int>::iterator
             coflow_in_q = current_queue.begin();
         coflow_in_q != current_queue.end();
      /*advance in loop */) {
      if (!ContainsKey(active_coflow_id_ptr_map, *coflow_in_q)) {
        // remove obsolete coflow
        coflow_in_q = current_queue.erase(coflow_in_q);
        continue;
      }
      // check if the coflow exceed the queue height.
      int coflow_jid = *coflow_in_q;
      Coflow *cfp = active_coflow_id_ptr_map[coflow_jid];
      double coflow_sent_byte = cfp->GetSentByte();
      unsigned int q_to_go = 0;
      for (double cutoff = AALO_INIT_Q_HEIGHT;
           cutoff < coflow_sent_byte;
           cutoff *= AALO_Q_HEIGHT_MULTI) {
        q_to_go++;
      }
      if (q_to_go != queue_runner) {
        if (DEBUG_LEVEL >= 1) {
          cout << "Qmv Q-" << queue_runner << " -> Q-" << q_to_go
               << " bytes " << cfp->GetSentByte() << "/" << cfp->GetSizeInByte()
               << " " << cfp->toString() << endl;
        }
        // add to target queue
        coflow_id_queues[q_to_go].push_back(coflow_jid);
        // rm from current queue
        coflow_in_q = current_queue.erase(coflow_in_q);
        continue;
      }

      // otherwise, no updateing to current coflow
      // so we advance the iterator
      coflow_in_q++;
      continue;
    }
  }

  if (DEBUG_LEVEL >= 5) {
    for (unsigned int queue_runner = 0;
         queue_runner < AALO_Q_NUM;
         queue_runner++) {
      vector<int> &current_queue = coflow_id_queues[queue_runner];
      cout << "Q-" << queue_runner << " : ";
      for (vector<int>::iterator
               coflow_in_q = current_queue.begin();
           coflow_in_q != current_queue.end();
           coflow_in_q++) {
        cout << *coflow_in_q << " ";
      }
      cout << endl;
    }
  }

}