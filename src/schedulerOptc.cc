//
//  schedulerOptc.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 2/5/15.
//  Copyright (c) 2015 Xin Sunny Huang. All rights reserved.
//

#include <algorithm>
#include <iomanip>
#include <sys/time.h>
#include <sstream>
#include <string.h>

#include "scheduler.h"
#include "events.h"
#include "global.h"
#include "util.h"
#include "coflow.h"
#include "hopcroft.h"
#include "edmond.h"

#define MWM_RANGE MAXWT //2147483647 = 2,147,483,647

///////////////////////////////////////////////////////
////////////// Code for base class - Opct scheduler
///////////////////////////////////////////////////////

bool slotLengthComp(SLOT *s1, SLOT *s2) {
  return s1->timeLength > s2->timeLength;
}

ostream &operator<<(std::ostream &out, const PortStatus &value) {
  switch (value) {
    case kPortUnknown:return out << "Unknown";
    case kPortTX:return out << "TX";
    case kPortIdle:return out << "Idle";
    case kPortSolidSwitch:return out << "SolidSwitch";
    case kPortFakeSwitch:return out << "FakeSwitch";
    case kPortWait:return out << "Wait";
    default:break;
  }
  return out << "";
}

SchedulerOptc::SchedulerOptc() : Scheduler() {
  m_max_optimal_workspan_src = -1;
  m_max_optimal_workspan_dst = -1;

  m_circuitAuditFile.open(CIRCUIT_AUDIT_FILE_NAME);
  if (!m_circuitAuditFile.is_open()) {
    cout << "Error: unable to open circuit audit file "
         << CIRCUIT_AUDIT_FILE_NAME << endl;
    cout << "Now terminate the program" << endl;
    exit(-1);
  }
}

SchedulerOptc::~SchedulerOptc() {
  if (m_circuitAuditFile.is_open()) {
    m_circuitAuditFile.close();
  }
}

void
SchedulerOptc::GetDemandAll(map<pair<int, int>, long> &demandMatrix) {
  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {
    // for each coflow

    if ((*cfIt)->IsComplete()) {
      //such coflow has completed
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
      MapWithInc(demandMatrix,
                 make_pair((*fpIt)->GetSrc(), (*fpIt)->GetDest()),
                 (*fpIt)->GetBitsLeft());
    }
  }
}

void
SchedulerOptc::DemandToDoublySM(map<pair<int, int>, long> &demandMatrix) {

  if (!FILL_DEMAND) return;

  map<int, long> sBpsUsed;
  map<int, long> rBpsUsed;
  map<int, long> sBitFree;
  map<int, long> rBitFree;


  // STEP 1: calculate sBitFree, rBitFree, bitBudget;
  for (map<pair<int, int>, long>::iterator dmIt = demandMatrix.begin();
       dmIt != demandMatrix.end(); dmIt++) {
    MapWithInc(sBpsUsed, dmIt->first.first, dmIt->second);
    MapWithInc(rBpsUsed, dmIt->first.second, dmIt->second);
  }

  long bitBudget = 0;
  bitBudget = max(bitBudget, MaxMap(sBpsUsed));
  bitBudget = max(bitBudget, MaxMap(rBpsUsed));

  for (map<int, long>::iterator sUsedIt = sBpsUsed.begin();
       sUsedIt != sBpsUsed.end(); sUsedIt++) {
    MapWithDef(sBitFree, sUsedIt->first, bitBudget - sUsedIt->second);
  }

  for (map<int, long>::iterator rUsedIt = rBpsUsed.begin();
       rUsedIt != rBpsUsed.end(); rUsedIt++) {
    MapWithDef(rBitFree, rUsedIt->first, bitBudget - rUsedIt->second);
  }

  long numSrc = sBitFree.size();
  long numDst = rBitFree.size();
  long numEnd = max(numSrc, numDst);

  // Step 2: fill with dummy src/dsr for a SQUARE matrix
  // filling is only reflected in sBitFree/rBitFree
  // demand matrix unmodified.
  // 1/2 either one or none of the branch to be exe to fill demand
  if (numEnd > numSrc) {

    // fill out with (numEnd - numSrc) randomly selected src
    long numSrc2Add = numEnd - numSrc;

    while (numSrc2Add > 0) {
      int randSrc = (rand() % NUM_RACK);
      if (0 > FindWithDef(sBitFree, randSrc, (long) -1)) {
        // src available
        MapWithDef(sBitFree, randSrc, bitBudget);
        //added successfully
        numSrc2Add--;
        //debug
        // cout << "add src " << randSrc << endl;
      }
    }
  }
  // 2/2 either one or none of the branch to be exe to fill demand
  if (numEnd > numDst) {

    // fill out with (numEnd - numSrc) randomly selected dst
    long numDst2Add = numEnd - numDst;

    while (numDst2Add > 0) {
      int randDst = (rand() % NUM_RACK);
      if (0 > FindWithDef(rBitFree, randDst, (long) -1)) {
        // dst available
        MapWithDef(rBitFree, randDst, bitBudget);
        //added successfully
        numDst2Add--;
        //debug
        // cout << "add dst " << randDst << endl;
      }
    }
  }

  // STEP 2: call another DemandToDoublySM with
  //                  sBitFree, rBitFree, bitBudget;
  // DemandToDoublySM(demandMatrix, sBitFree, rBitFree, bitBudget);
  SolsticeQuickStuff(demandMatrix, sBitFree, rBitFree, bitBudget);
}

// trying to modified demand matrix.
// suboptimal: potentially making the matrix less sparse
void
SchedulerOptc::FillSquareMatrixToDSM(map<pair<int, int>, long> &demandMatrix,
                                     map<int, long> &sBitFree,
                                     map<int, long> &rBitFree,
                                     long &bitBudget) {

  if (!FILL_DEMAND) return;

  map<int, long> sBitFree2Fill = sBitFree;
  map<int, long> rBitFree2Fill = rBitFree;

  // now that we have a SQUARE matrix
  // fill demand into a SQUARE full DOUBLY stochastic matrix

  while (1) {
    // TODO: there might be other way to pick up a victim
    // here we simply pick the src and dst that has
    // the max free resource
    // or in other word, we pick the src and dst that has
    // the least demand
    int minSrc = Key2MaxPositiveMap(sBitFree2Fill);
    int minDst = Key2MaxPositiveMap(rBitFree2Fill);
    if (sBitFree2Fill[minSrc] <= 0 && rBitFree2Fill[minDst] <= 0) {
      break;
    } else if ((sBitFree2Fill[minSrc] <= 0) != (rBitFree2Fill[minDst] <= 0)) {
      cout << "demand matrix filling error : "
           << "sBitFree2Fill should reach zero at the same time" << endl;
      break;
    }
    //inc should always be larger than zero
    long inc = min(sBitFree2Fill[minSrc], rBitFree2Fill[minDst]);
    MapWithInc(demandMatrix,
               make_pair(minSrc, minDst),
               inc);
    sBitFree2Fill[minSrc] -= inc;
    rBitFree2Fill[minDst] -= inc;
    //cout <<  "fill " << minSrc << "->" << minDst << " with " << inc<< endl;
  }

  //valid
  ValidateDSM(demandMatrix);
}

// trying to modified demand matrix.
// a better way to modify demand matrix by maintaining sparsity.
void
SchedulerOptc::SolsticeQuickStuff(map<pair<int, int>, long> &demandMatrix,
                                  map<int, long> &sBitFree,
                                  map<int, long> &rBitFree,
                                  long &bitBudget) {
  if (!FILL_DEMAND) return;

  // now that we have a SQUARE matrix
  // fill demand into a SQUARE full DOUBLY stochastic matrix

  map<int, long> sBitFree2Fill = sBitFree;
  map<int, long> rBitFree2Fill = rBitFree;

  // for non-zero entries
  for (map<pair<int, int>, long>::iterator
           demand_iter = demandMatrix.begin();
       demand_iter != demandMatrix.end();
       demand_iter++) {
    if (demand_iter->second <= 0) {
      continue;
    }
    int src = demand_iter->first.first;
    int dst = demand_iter->first.second;
    long src_free_bit = MapWithDef(sBitFree2Fill, src, (long) 0);
    long dst_free_bit = MapWithDef(rBitFree2Fill, dst, (long) 0);
    long bit_to_add = min(src_free_bit, dst_free_bit);
    if (bit_to_add > 0) {
      demand_iter->second += bit_to_add;
      sBitFree2Fill[src] -= bit_to_add;
      rBitFree2Fill[dst] -= bit_to_add;
      if (DEBUG_LEVEL >= 15) {
        cout << "Non-zero Stuffing " << src << "->" << dst
             << ":" << bit_to_add << " " << bit_to_add / (double) OPTC_BPS
             << endl;
      }
    }
  }

  // for zero entries.
  for (map<int, long>::iterator
           src_pair = sBitFree2Fill.begin();
       src_pair != sBitFree2Fill.end();
       src_pair++) {
    int src = src_pair->first;
    for (map<int, long>::iterator
             dst_pair = rBitFree2Fill.begin();
         dst_pair != rBitFree2Fill.end();
         dst_pair++) {
      int dst = dst_pair->first;
      pair<int, int> src_dst_pair(src, dst);
      if (!ContainsKey(demandMatrix, src_dst_pair)) {
        // this is a zero entry
        long src_free_bit = MapWithDef(sBitFree2Fill, src, (long) 0);
        long dst_free_bit = MapWithDef(rBitFree2Fill, dst, (long) 0);
        long bit_to_add = min(src_free_bit, dst_free_bit);
        if (bit_to_add > 0) {
          demandMatrix.insert(std::make_pair(src_dst_pair, bit_to_add));
          sBitFree2Fill[src] -= bit_to_add;
          rBitFree2Fill[dst] -= bit_to_add;
          if (DEBUG_LEVEL >= 15) {
            cout << "Zero Stuffing " << src << "->" << dst
                 << ":" << bit_to_add << " " << bit_to_add / (double) OPTC_BPS
                 << endl;
          }
        }

      }
    }
  }
  //valid
  ValidateDSM(demandMatrix);
}

bool
SchedulerOptc::ValidateCircuitAssignment(const set<pair<int, int>> &circuits) {
  set<int> busy_src;
  set<int> busy_dst;
  for (set<pair<int, int>>::const_iterator circuit = circuits.begin();
       circuit != circuits.end(); circuit++) {
    if (ContainsKey(busy_src, circuit->first)
        || ContainsKey(busy_dst, circuit->second)) {
      cout << "circuit assignment error: "
           << " src port " << circuit->first << " and/or "
           << " dst port " << circuit->second
           << " has been assigned!" << endl;
      return false;
    } else {
      busy_src.insert(circuit->first);
      busy_dst.insert(circuit->second);
    }
  }
  // cout << "Circuit assignment is valid." << endl;
  return true;
}

void
SchedulerOptc::ValidateDSM(map<pair<int, int>, long> &demandMatrix) {

  if (!VALIDATE_DSM) return;

  if (DEBUG_LEVEL >= 5) {
    cout << " demand matrix " << endl;
    for (map<pair<int, int>, long>::const_iterator
             dm_iter = demandMatrix.begin();
         dm_iter != demandMatrix.end();
         dm_iter++) {
      cout << dm_iter->first.first << "->" << dm_iter->first.second
           << ":" << dm_iter->second << endl;
    }
  }

  map<int, long> sValidate;
  map<int, long> rValidate;

  for (map<pair<int, int>, long>::iterator dmIt = demandMatrix.begin();
       dmIt != demandMatrix.end(); dmIt++) {
    int src = dmIt->first.first;
    int dst = dmIt->first.second;
    long dm = dmIt->second;
    MapWithInc(sValidate, src, dm);
    MapWithInc(rValidate, dst, dm);
  }
  if (sValidate.size() != rValidate.size()) {
    cout << "DSM error: matrix after fill is not square " << endl;
    cout << "src num = " << sValidate.size() << endl;
    cout << "dst num = " << rValidate.size() << endl;
    return;
  }
  if (MaxMap(sValidate) != MaxMap(rValidate)) {
    cout << "DSM error: row sum and col sum not equal " << endl;
    cout << "max src sum = " << MaxMap(sValidate) << endl;
    cout << "max dst sum = " << MaxMap(rValidate) << endl;
    return;
  }

  long bitBudget = MaxMap(sValidate);
  for (map<int, long>::iterator srcIt = sValidate.begin();
       srcIt != sValidate.end(); srcIt++) {
    if (srcIt->second != bitBudget) {
      cout << "DSM error: src " << srcIt->first
           << " sum " << srcIt->second
           << " not equal to bitBudget " << bitBudget << endl;
      return;
    }
  }
  for (map<int, long>::iterator dstIt = rValidate.begin();
       dstIt != rValidate.end(); dstIt++) {
    if (dstIt->second != bitBudget) {
      cout << "DSM error: dst " << dstIt->first
           << " sum " << dstIt->second
           << " not equal to bitBudget " << bitBudget << endl;
      return;
    }
  }

  cout << "DSM validated " << endl;
}

// mark down circuit activities
// in src/dst/flow activity.
// only does auditing for the most prioritized coflow.
// bottleneck src/dst, demanding src/dst/circuit
// ONLY apply to the first coflow.
// ideally, works under coflow 1by1 traffic generator.
// required: only 1 coflow running.
void
SchedulerOptc::CircuitAuditIfNeeded(double beginTime,
                                    double endTime) {

  if (!PORT_AUDIT_ALL_SRC
      && !PORT_AUDIT_BOTTLENECK
      && !PORT_AUDIT_LAST_GONE_SRC_DST) {
    // return directly if not doing any port auditing.
    return;
  }

  if (beginTime >= endTime) {
    return;
  }

  vector<map<pair<int, int>, long>> layered_demand;
  CoflowsToLayeredDemands(m_coflowPtrVector, layered_demand);
  if (layered_demand.empty()) {
    return;
  }
  map<pair<int, int>, long> first_coflow_demand = layered_demand[0];

  // find out bottleneck src/dst if needed.
  if (m_max_optimal_workspan_src <= -1 && m_max_optimal_workspan_dst <= -1) {
    if (m_coflowPtrVector.size() != 1) {
      cout
          << " Error when trying to find out the first coflow on the coflow list."
          << endl;
    } else {
      // safe to access
      // find out the port with max optimal-work-span
      // (with switching and payload considered)
      m_coflowPtrVector[0]->GetPortOnMaxOptimalWorkSpan(
          m_max_optimal_workspan_src,
          m_max_optimal_workspan_dst);
      if (DEBUG_LEVEL >= 3) {
        cout << "max optimal work span on (only one of src/dst valid)"
             << " src " << m_max_optimal_workspan_src
             << " dst " << m_max_optimal_workspan_dst << endl;
      }
    }
  }

  // find out src/dst with non-zero demand.
  set<int> demanding_src;
  set<int> demanding_dst;
  for (map<pair<int, int>, long>::const_iterator
           demand_pair = first_coflow_demand.begin();
       demand_pair != first_coflow_demand.end();
       demand_pair++) {
    demanding_src.insert(demand_pair->first.first);
    demanding_dst.insert(demand_pair->first.second);
  }

  // for src/dst that is currently working.
  set<pair<int, int>> active_circuits = m_activeCircuits;
  set<int> active_src;
  set<int> active_dst;
  map<int, int> active_src_to_dst;
  map<int, int> active_dst_to_src;
  for (set<pair<int, int>>::const_iterator
           circuit = active_circuits.begin();
       circuit != active_circuits.end();
       circuit++) {
    // active_src/active_dst logics can be further collapesed
    // by using active_src_to_dst/active_dst_to_src only.
    active_src.insert(circuit->first);
    active_dst.insert(circuit->second);
    active_src_to_dst[circuit->first] = circuit->second;
    active_dst_to_src[circuit->second] = circuit->first;
  }

  // for src/dst that is currently switching.
  set<pair<int, int>> switching_circuits = m_circuits_switching_to;
  set<int> switching_src;
  set<int> switching_dst;
  map<int, int> switching_src_to_dst;
  map<int, int> switching_dst_to_src;
  for (set<pair<int, int>>::const_iterator
           circuit = switching_circuits.begin();
       circuit != switching_circuits.end();
       circuit++) {
    switching_src.insert(circuit->first);
    switching_dst.insert(circuit->second);
    switching_src_to_dst[circuit->first] = circuit->second;
    switching_dst_to_src[circuit->second] = circuit->first;
  }

  // log each demanding src.
  for (set<int>::const_iterator
           src_iter = demanding_src.begin();
       src_iter != demanding_src.end();
       src_iter++) {
    int src = *src_iter;
    pair<int, int> circuit_if_active;
    if (ContainsKey(active_src, src)) {
      // port up.
      circuit_if_active = pair<int, int>(src, active_src_to_dst[src]);
    } else {
      circuit_if_active = pair<int, int>(-1, -1);
    }
    pair<int, int> circuit_if_switching;
    if (ContainsKey(switching_src, src)) {
      // port switching.
      circuit_if_switching = pair<int, int>(src, switching_src_to_dst[src]);
    } else {
      circuit_if_switching = pair<int, int>(-1, -1);
    }
    PortStatus src_status = CheckPortStatus(src, true,
                                            circuit_if_active,
                                            circuit_if_switching,
                                            first_coflow_demand);
    if (src_status == kPortUnknown) {
      cout << " Warming: see unknown port status for src " << src << "!"
           << endl;
    }
    m_src_activity[src].push(PortActivity(beginTime,
                                          endTime,
                                          src_status));
    if (DEBUG_LEVEL >= 9) {
      cout << "src " << src << " " << src_status
           << " (" << beginTime << ", " << endTime << ")" << endl;
    }
  }

  // log each demanding dst.
  for (set<int>::const_iterator
           dst_iter = demanding_dst.begin();
       dst_iter != demanding_dst.end();
       dst_iter++) {
    int dst = *dst_iter;
    pair<int, int> circuit_if_active;
    if (ContainsKey(active_dst, dst)) {
      // port up.
      circuit_if_active = pair<int, int>(active_dst_to_src[dst], dst);
    } else {
      circuit_if_active = pair<int, int>(-1, -1);
    }
    pair<int, int> circuit_if_switching;
    if (ContainsKey(switching_dst, dst)) {
      // port switching.
      circuit_if_switching = pair<int, int>(switching_dst_to_src[dst], dst);
    } else {
      circuit_if_switching = pair<int, int>(-1, -1);
    }
    PortStatus dst_status = CheckPortStatus(dst, false,
                                            circuit_if_active,
                                            circuit_if_switching,
                                            first_coflow_demand);
    if (dst_status == kPortUnknown) {
      cout << " Warming: see unknown port status for dst " << dst << "!"
           << endl;
    }
    m_dst_activity[dst].push(PortActivity(beginTime,
                                          endTime,
                                          dst_status));
    if (DEBUG_LEVEL >= 9) {
      cout << "dst " << dst << " " << dst_status
           << " (" << beginTime << ", " << endTime << ")" << endl;
    }
  }

  // log bottleneck port.
  {
    bool bottleneck_on_src =
        m_max_optimal_workspan_src > m_max_optimal_workspan_dst;
    int bottleneck_port = bottleneck_on_src ? m_max_optimal_workspan_src
                                            : m_max_optimal_workspan_dst;
    pair<int, int> circuit_if_active;
    if (bottleneck_on_src && ContainsKey(active_src, bottleneck_port)) {
      // bottleneck on src.
      int src = bottleneck_port;
      circuit_if_active = pair<int, int>(src, active_src_to_dst[src]);
    } else if (!bottleneck_on_src
        && ContainsKey(active_dst, bottleneck_port)) {
      // bottleneck on dst.
      int dst = bottleneck_port;
      circuit_if_active = pair<int, int>(active_dst_to_src[dst], dst);
    } else {
      circuit_if_active = pair<int, int>(-1, -1);
    }
    pair<int, int> circuit_if_switching;
    if (bottleneck_on_src && ContainsKey(switching_src, bottleneck_port)) {
      // src port switching.
      int src = bottleneck_port;
      circuit_if_switching = pair<int, int>(src, switching_src_to_dst[src]);
    } else if (!bottleneck_on_src
        && ContainsKey(switching_dst, bottleneck_port)) {
      // dst port switching.
      int dst = bottleneck_port;
      circuit_if_switching = pair<int, int>(switching_dst_to_src[dst], dst);
    } else {
      circuit_if_switching = pair<int, int>(-1, -1);
    }
    PortStatus port_status = CheckPortStatus(bottleneck_port,
                                             bottleneck_on_src,
                                             circuit_if_active,
                                             circuit_if_switching,
                                             first_coflow_demand);
    if (port_status == kPortUnknown) {
      string port_type = bottleneck_on_src ? " src " : " dst ";
      cout << " Warming: see unknown port status for " << port_type
           << bottleneck_port << "!" << endl;
    }
    m_max_optimal_workspan_port_activity.push(PortActivity(beginTime,
                                                           endTime,
                                                           port_status));
    if (DEBUG_LEVEL >= 9) {
      string port_type = bottleneck_on_src ? " src " : " dst ";
      cout << "bottleneck " << port_type << bottleneck_port << " "
           << port_status << " (" << beginTime << ", " << endTime << ")"
           << endl;
    }
  }

}

// when port is active, both entries on circuit_if_active >= 0
// if either entry on circuit_if_active <0 : port not active
// when port is switching, both entries on circuit_if_switching >= 0
// if either entry on circuit_if_switching <0 : port not switching
PortStatus
SchedulerOptc::CheckPortStatus(int port, bool is_src,
                               const pair<int, int> &circuit_if_active,
                               const pair<int, int> &circuit_if_switching,
                               const map<pair<int, int>, long> &demand) {

  bool is_on_active_circuit =
      circuit_if_active.first >= 0 && circuit_if_active.second >= 0;
  bool is_switching =
      circuit_if_switching.first >= 0 && circuit_if_switching.second >= 0;
  if (is_on_active_circuit) {
    // port up.
    if (ContainsKey(demand, circuit_if_active)) {
      // port up and serving demand
      return kPortTX;
    } else {
      // port up but not serving demand
      return kPortIdle;
    }
  } else {
    // port down.
    if (is_switching) {
      // port down due to switching.
      if (ContainsKey(demand, circuit_if_switching)) {
        // switching for something.
        return kPortSolidSwitch;
      } else {
        // switching for nothing: fake.
        return kPortFakeSwitch;
      }
    } else {
      // port down but not in switching => waiting.
      return kPortWait;
    }
  }

  return kPortUnknown;
}

// do some calculation and
// write out circuit statistics.
void
SchedulerOptc::WriteCircuitAuditIfNeeded(double endTime,
                                         vector<Coflow *> &coflows_done) {

  if (!PORT_AUDIT_ALL_SRC
      && !PORT_AUDIT_BOTTLENECK
      && !PORT_AUDIT_LAST_GONE_SRC_DST) {
    // return directly if not doing any port auditing.
    return;
  }

  if (coflows_done.size() != 1) {
    cout << "Error in port auditing: expecting a coflow to finish only."
         << endl;
    return;
  }

  // calculate on coflow
  map<pair<int, int>, long> first_coflow_demand_bit;
  vector<Flow *> *flow_vec_ptr = coflows_done[0]->GetFlows();
  int job_id = coflows_done[0]->GetJobId();
  for (vector<Flow *>::const_iterator
           fpIt = flow_vec_ptr->begin();
       fpIt != flow_vec_ptr->end();
       fpIt++) {
    pair<int, int> circuit_to_take((*fpIt)->GetSrc(), (*fpIt)->GetDest());
    first_coflow_demand_bit.insert(std::make_pair(circuit_to_take,
                                                  (*fpIt)->GetSizeInBit()));
  }
  double this_coflow_CCT = endTime - coflows_done[0]->GetStartTime();
  double this_coflow_bottleneck_sec =
      coflows_done[0]->GetMaxOptimalWorkSpanInSeconds();
  double this_coflow_bottleneck_payload_sec =
      coflows_done[0]->GetLoadOnMaxOptimalWorkSpanInBits() / (double) OPTC_BPS;

  // calculate on last_gone src/dst.
  if (PORT_AUDIT_LAST_GONE_SRC_DST) {
    // find last_gone src.
    double last_src_gone_time = -1;
    int last_src_gone = -1;
    for (PortActivityMap::const_iterator
             src_act_iter = m_src_activity.begin();
         src_act_iter != m_src_activity.end();
         src_act_iter++) {
      PortActivity src_last_act = src_act_iter->second.top();
      int src = src_act_iter->first;
      if (last_src_gone_time < src_last_act.m_end_time
          || last_src_gone_time < 0) {
        last_src_gone_time = src_last_act.m_end_time;
        last_src_gone = src;
      }
    }

    // find last_gone dst.
    double last_dst_gone_time = -1;
    int last_dst_gone = -1;
    for (PortActivityMap::const_iterator
             dst_act_iter = m_dst_activity.begin();
         dst_act_iter != m_dst_activity.end();
         dst_act_iter++) {
      PortActivity dst_last_act = dst_act_iter->second.top();
      int dst = dst_act_iter->first;
      if (last_dst_gone_time < dst_last_act.m_end_time
          || last_dst_gone_time < 0) {
        last_dst_gone_time = dst_last_act.m_end_time;
        last_dst_gone = dst;
      }
    }

    PortActivityQueue last_gone_src_act_queue = m_src_activity[last_src_gone];
    PortActivityQueue last_gone_dst_act_queue = m_dst_activity[last_dst_gone];

    // calculate last_gone src activity stats.
    double last_gone_src_solid_switch_time_sum =
        0.0;  // down for switching and switching for demand
    double last_gone_src_tx_time_sum = 0.0;            // up and tx
    double last_gone_src_wait_time_sum =
        0.0;          // wasted: just hang there. not in active circuit. not switching.
    double last_gone_src_fake_switch_time_sum =
        0.0;   // wasted: down for switching but switching for no demand
    double last_gone_src_idle_time_sum = 0.0;          // wasted: up not tx
    while (!last_gone_src_act_queue.empty()) {
      PortActivity src_act = last_gone_src_act_queue.top();
      if (kPortTX == src_act.m_status) {
        last_gone_src_tx_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      } else if (kPortSolidSwitch == src_act.m_status) {
        last_gone_src_solid_switch_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      } else if (kPortWait == src_act.m_status) {
        last_gone_src_wait_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      } else if (kPortFakeSwitch == src_act.m_status) {
        last_gone_src_fake_switch_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      } else if (kPortIdle == src_act.m_status) {
        last_gone_src_idle_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      }
      last_gone_src_act_queue.pop();
    }

    double last_gone_dst_solid_switch_time_sum =
        0.0;  // down for switching and switching for demand
    double last_gone_dst_tx_time_sum = 0.0;            // up and tx
    double last_gone_dst_wait_time_sum =
        0.0;          // wasted: just hang there. not in active circuit. not switching.
    double last_gone_dst_fake_switch_time_sum =
        0.0;   // wasted: down for switching but switching for no demand
    double last_gone_dst_idle_time_sum = 0.0;          // wasted: up not tx
    while (!last_gone_dst_act_queue.empty()) {
      PortActivity dst_act = last_gone_dst_act_queue.top();
      if (kPortTX == dst_act.m_status) {
        last_gone_dst_tx_time_sum +=
            (dst_act.m_end_time - dst_act.m_start_time);
      } else if (kPortSolidSwitch == dst_act.m_status) {
        last_gone_dst_solid_switch_time_sum +=
            (dst_act.m_end_time - dst_act.m_start_time);
      } else if (kPortWait == dst_act.m_status) {
        last_gone_dst_wait_time_sum +=
            (dst_act.m_end_time - dst_act.m_start_time);
      } else if (kPortFakeSwitch == dst_act.m_status) {
        last_gone_dst_fake_switch_time_sum +=
            (dst_act.m_end_time - dst_act.m_start_time);
      } else if (kPortIdle == dst_act.m_status) {
        last_gone_dst_idle_time_sum +=
            (dst_act.m_end_time - dst_act.m_start_time);
      }
      last_gone_dst_act_queue.pop();
    }

    // calculate last_gone src/dst payload.
    double last_gone_src_payload_sec
        = coflows_done[0]->GetLoadOnPortInBits(last_src_gone, -1)
            / (double) OPTC_BPS;
    double last_gone_dst_payload_sec
        = coflows_done[0]->GetLoadOnPortInBits(-1, last_dst_gone)
            / (double) OPTC_BPS;



    // write stats on last_gone src/dst.
    if (!m_port_audit_title_line_last_gone) {
      m_port_audit_title_line_last_gone = true;

      m_portAuditFile
          << "jobid" << '\t'
          << "cct" << '\t'
          << "bottleload" << '\t'
          << "bottletime" << '\t'
          << "src_load" << '\t'
          << "src_tx" << '\t'
          << "src_solid_switch" << '\t'
          << "src_fake_switch" << '\t'
          << "src_idle" << '\t'
          << "src_wait" << '\t'
          << "dst_load" << '\t'
          << "dst_tx" << '\t'
          << "dst_solid_switch" << '\t'
          << "dst_fake_switch" << '\t'
          << "dst_idle" << '\t'
          << "dst_wait" << '\t'
          << endl;
    }
    m_portAuditFile
        << job_id << '\t'
        << this_coflow_CCT << '\t'
        << this_coflow_bottleneck_payload_sec << '\t'
        << this_coflow_bottleneck_sec << '\t'
        << last_gone_src_payload_sec << '\t' // src
        << last_gone_src_tx_time_sum << '\t'
        << last_gone_src_solid_switch_time_sum << '\t'
        << last_gone_src_fake_switch_time_sum << '\t'
        << last_gone_src_idle_time_sum << '\t'
        << last_gone_src_wait_time_sum << '\t'
        << last_gone_dst_payload_sec << '\t' // dst
        << last_gone_dst_tx_time_sum << '\t'
        << last_gone_dst_solid_switch_time_sum << '\t'
        << last_gone_dst_fake_switch_time_sum << '\t'
        << last_gone_dst_idle_time_sum << '\t'
        << last_gone_dst_wait_time_sum << '\t'
        << endl;

  } // if (PORT_AUDIT_LAST_GONE_SRC_DST)

  // calculate on bottleneck src/dst.
  if (PORT_AUDIT_BOTTLENECK) {
    bool bottleneck_on_src =
        m_max_optimal_workspan_src > m_max_optimal_workspan_dst;
    int bottleneck_port = bottleneck_on_src ? m_max_optimal_workspan_src
                                            : m_max_optimal_workspan_dst;
    string bottleneck_port_type = bottleneck_on_src ? " src " : " dst ";
    cout << "working on bottleneck " << bottleneck_port_type << bottleneck_port
         << endl;

    PortActivityQueue
        bottleneck_port_act_queue = m_max_optimal_workspan_port_activity;
    double bottleneck_port_solid_switch_time_sum =
        0.0;  // down for switching and switching for demand
    double bottleneck_port_tx_time_sum = 0.0;            // up and tx
    double
        bottleneck_port_wait_time_sum =
        0.0;          // wasted: just hang there. not in active circuit. not switching.
    double bottleneck_port_fake_switch_time_sum =
        0.0;   // wasted: down for switching but switching for no demand
    double bottleneck_port_idle_time_sum = 0.0;          // wasted: up not tx
    while (!bottleneck_port_act_queue.empty()) {
      PortActivity src_act = bottleneck_port_act_queue.top();
      if (kPortTX == src_act.m_status) {
        bottleneck_port_tx_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      } else if (kPortSolidSwitch == src_act.m_status) {
        bottleneck_port_solid_switch_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      } else if (kPortWait == src_act.m_status) {
        bottleneck_port_wait_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      } else if (kPortFakeSwitch == src_act.m_status) {
        bottleneck_port_fake_switch_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      } else if (kPortIdle == src_act.m_status) {
        bottleneck_port_idle_time_sum +=
            (src_act.m_end_time - src_act.m_start_time);
      }
      bottleneck_port_act_queue.pop();
    }

    // calculate bottleneck src/dst payload.
    double bottleneck_port_payload_sec =
        coflows_done[0]->GetLoadOnPortInBits(m_max_optimal_workspan_src,
                                             m_max_optimal_workspan_dst)
            / (double) OPTC_BPS;

    // write stats on with max-optimal-work-span src/dst.
    if (!m_port_audit_title_line_bottleneck) {
      m_port_audit_title_line_bottleneck = true;
      m_portAuditFile
          << "jobid" << '\t'
          << "cct" << '\t'
          << "bottleload" << '\t'
          << "bottletime" << '\t'
          << "bottle_load" << '\t'
          << "bottle_tx" << '\t'
          << "bottle_solid_switch" << '\t'
          << "bottle_fake_switch" << '\t'
          << "bottle_idle" << '\t'
          << "bottle_wait" << '\t'
          << endl;
    }
    m_portAuditFile << setw(8) << job_id << '\t'
                    << this_coflow_CCT << '\t'
                    << this_coflow_bottleneck_payload_sec << '\t'
                    << this_coflow_bottleneck_sec << '\t'
                    << bottleneck_port_payload_sec << '\t'
                    << bottleneck_port_tx_time_sum << '\t'
                    << bottleneck_port_solid_switch_time_sum << '\t'
                    << bottleneck_port_fake_switch_time_sum << '\t'
                    << bottleneck_port_idle_time_sum << '\t'
                    << bottleneck_port_wait_time_sum << '\t'
                    << endl;

  } // if (PORT_AUDIT_BOTTLENECK)

  // calculate on the whole exe graph.
  // based on src port.
  if (PORT_AUDIT_ALL_SRC) {
    double all_src_solid_switch_time_sum =
        0.0;  // down for switching and switching for demand
    double all_src_tx_time_sum = 0.0;            // up and tx
    double all_src_wait_time_sum =
        0.0;          // wasted: just hang there. not in active circuit. not switching.
    double all_src_fake_switch_time_sum =
        0.0;   // wasted: down for switching but switching for no demand
    double all_src_idle_time_sum = 0.0;          // wasted: up not tx

    for (PortActivityMap::const_iterator
             src_act_iter = m_src_activity.begin();
         src_act_iter != m_src_activity.end();
         src_act_iter++) {
      PortActivityQueue src_act_queue = src_act_iter->second;
      while (!src_act_queue.empty()) {
        PortActivity src_act = src_act_queue.top();
        if (kPortTX == src_act.m_status) {
          all_src_tx_time_sum += (src_act.m_end_time - src_act.m_start_time);
        } else if (kPortSolidSwitch == src_act.m_status) {
          all_src_solid_switch_time_sum +=
              (src_act.m_end_time - src_act.m_start_time);
        } else if (kPortWait == src_act.m_status) {
          all_src_wait_time_sum += (src_act.m_end_time - src_act.m_start_time);
        } else if (kPortFakeSwitch == src_act.m_status) {
          all_src_fake_switch_time_sum +=
              (src_act.m_end_time - src_act.m_start_time);
        } else if (kPortIdle == src_act.m_status) {
          all_src_idle_time_sum += (src_act.m_end_time - src_act.m_start_time);
        }
        src_act_queue.pop();
      }
    }
    unsigned long this_coflow_flow_num = coflows_done[0]->GetFlows()->size();
    double all_src_payload_sec =
        coflows_done[0]->GetSizeInByte() * 8.0 / (double) OPTC_BPS;
    if (!m_port_audit_title_line_all) {
      m_port_audit_title_line_all = true;
      m_portAuditFile
          << "jobid" << '\t'
          << "cct" << '\t'
          << "flow#" << '\t'
          << "all_solid_switch" << '\t'
          << "bottletime" << '\t'
          << "all_load" << '\t'
          << "all_tx" << '\t'
          << "all_solid_switch" << '\t'
          << "all_fake_switch" << '\t'
          << "all_idle" << '\t'
          << "all_wait" << '\t'
          << endl;
    }
    m_portAuditFile
        << job_id << '\t'
        << this_coflow_CCT << '\t'
        << this_coflow_flow_num << '\t'
        << all_src_solid_switch_time_sum << '\t'
        << this_coflow_bottleneck_sec << '\t'
        << all_src_payload_sec << '\t'
        << all_src_tx_time_sum << '\t'
        << all_src_solid_switch_time_sum << '\t'
        << all_src_fake_switch_time_sum << '\t'
        << all_src_idle_time_sum << '\t'
        << all_src_wait_time_sum << '\t'
        << endl;
  } // if (PORT_AUDIT_ALL_SRC)

  // End of auditing.
  // clean up before we move on.
  m_src_activity.clear();
  m_dst_activity.clear();
  // reset bottleneck src/dst so that
  // auditor function knows to update
  // bottleneck src/dst when new coflow arrives.
  m_max_optimal_workspan_src = -1;
  m_max_optimal_workspan_dst = -1;
  m_max_optimal_workspan_port_activity = PortActivityQueue();
}

// Generate layered demand matrix.
// A demand matrix on a layer is the aggregated demand for a coflow.
void
SchedulerOptc::CoflowsToLayeredDemands(vector<Coflow *> &coflows,
                                       vector<map<pair<int, int>,
                                                  long>> &layered_demand) {
  vector<int> coflow_ids_to_be_dump; // useless.
  CoflowsToLayeredDemands(coflows, layered_demand, coflow_ids_to_be_dump);

}

// calculate layered demand
// so that a layer represent the demand from a coflow
// excluding circuits with more prioritized coflow presented (predecessors).
void
SchedulerOptc::CoflowsToLayeredDemands(vector<Coflow *> &coflows,
                                       vector<map<pair<int, int>,
                                                  long>> &layered_demand,
                                       vector<int> &layered_coflow_ids) {

  set<pair<int, int>> circuits_with_predecessor;

  for (vector<Coflow *>::iterator cfIt = coflows.begin();
       cfIt != coflows.end(); cfIt++) {
    // for each coflow

    if ((*cfIt)->IsComplete() || 0 == (*cfIt)->GetAlpha()) {
      // such coflow has completed
      // or has zero demand
      continue;
    }
    layered_demand.push_back(map<pair<int, int>, long>());
    layered_coflow_ids.push_back((*cfIt)->GetCoflowId());
    map<pair<int, int>, long> &one_coflow_demand = layered_demand.back();

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();
    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }
      pair<int, int> circuit_to_take((*fpIt)->GetSrc(), (*fpIt)->GetDest());
      // exclude circuits that have predecessors
      if (ContainsKey(circuits_with_predecessor, circuit_to_take)) {
        continue;
      }
      // yay! I am the first in this circuit.
      circuits_with_predecessor.insert(circuit_to_take);
      MapWithInc(one_coflow_demand,
                 circuit_to_take,
                 (*fpIt)->GetBitsLeft());

    }
  }
  // debug
  if (DEBUG_LEVEL >= 16) {
    if (layered_demand.size() != layered_coflow_ids.size()) {
      cout << "Error in " << __func__
           << " : layered_demand.size() != layered_coflow_ids.size()" << endl;
    }
    for (unsigned int coflow_runner = 0;
         coflow_runner < layered_demand.size();
         coflow_runner++) {
      map<pair<int, int>, long>
          &one_coflow_demand = layered_demand[coflow_runner];
      int coflow_id = layered_coflow_ids[coflow_runner];
      cout << "coflow [" << coflow_id << "] demand " << endl;
      for (map<pair<int, int>, long>::const_iterator
               demand = one_coflow_demand.begin();
           demand != one_coflow_demand.end();
           demand++) {
        cout << demand->first.first << " -> " << demand->first.second
             << " : " << demand->second << endl;
      }
    }
  }
}

// take in an arbitrary matrix
// and iteratively perform Hopcroft matching and min deduction
// NOTE: the matrix does not need to be SQUARE.
// when the matrix is a SQUARE doubly-stochastic matrix
// this is exactly BvN
void
SchedulerOptc::HopDecomp(map<pair<int, int>, long> &demandMatrix,
                         double timeSinceCalStart) {
  if (!m_nextCircuitVector.empty() && !demandMatrix.empty()) {
    cout << "try to generate new assignments "
         << " for non-empty next circuit vector!" << endl;
  }

  cout << "now generating patterns" << endl;

  map<pair<int, int>, long> origDemandMatrix = demandMatrix;

  struct timeval pattern_cal_start_time;
  struct timeval pattern_cal_end_time;

  while (demandMatrix.size() > 0) {

    set<pair<int, int>> *circuits = new set<pair<int, int>>();

    // STEP 1: get start time
    gettimeofday(&pattern_cal_start_time, NULL);

    // STEP 2: find a matching
    //******* Hopcroft matching begin ***********
    //                                          *
    DemandMatchingHop(demandMatrix, circuits);
    //                                          *
    //******* Hopcroft matching end *************

    gettimeofday(&pattern_cal_end_time, NULL);

    // STEP 3: find circuit assignment metrics

    // STEP 3.1: demand served
    //calculate the demand served by this pattern of circuits
    //note: use min demand among the matching
    long mMinDemand = -1;

    //******* find MIN entry begin ***************************
    //                                                       *
    for (set<pair<int, int>>::iterator setIt = circuits->begin();
         setIt != circuits->end(); setIt++) {
      int src = setIt->first;
      int dest = setIt->second;
      long demand = demandMatrix[make_pair(src, dest)];
      if (mMinDemand > demand || mMinDemand < 0) {
        mMinDemand = demand;
      }
    }
    //                                                       *
    //******* find MIN entry end ***************************



    // STEP 3.2: slot length based on demand served
    double slot_length = SecureFinishTime(mMinDemand, OPTC_BPS);

    // STEP 3.3: slot efficiency score
    double slot_score = GetSlotEfficiency(circuits, origDemandMatrix);

    // STEP 3.4: control plan delay
    double patternCalTime
        = secondPass(pattern_cal_end_time, pattern_cal_start_time);

    // test of pattern computation delay
    if (ZERO_COMP_TIME) {
      timeSinceCalStart = 0.0;
      patternCalTime = 0.0;
    }

    timeSinceCalStart += patternCalTime;
    double pattern_avail_time = m_currentTime + timeSinceCalStart;

    // STEP 3.5: a new assignment is generated!
    SLOT *slotp =
        new SLOT(slot_length, circuits, slot_score, pattern_avail_time);
    m_nextCircuitVector.push_back(slotp);

//        cout << " slot #" << m_nextCircuitVector.size()
//            << " length " << slotp->timeLength << endl;

    //Step 4: update demand matrix and remove zero demand
    for (set<pair<int, int>>::iterator setIt = circuits->begin();
         setIt != circuits->end(); setIt++) {
      int src = setIt->first;
      int dest = setIt->second;
      // use max demand
      demandMatrix[make_pair(src, dest)] -= mMinDemand;
      if (demandMatrix[make_pair(src, dest)] <= 0) {
        //remove entry
        // map<pair<int, int>, long>::iterator demand2RmIt;
        // demand2RmIt = demandMatrix.find(make_pair(src, dest));
        demandMatrix.erase(make_pair(src, dest));
      }
    }
  }
}

// takes a demand matrix as a graph
// calculate hopcroft (max cardinality) matching on the demand matrix
// stores the matching in the set designated by the circuit set pointer
// demand matrix should not be modified
void
SchedulerOptc::DemandMatchingHop(const map<pair<int, int>, long> &demandMatrix,
                                 set<pair<int, int>> *mSetPtr) {

  map<int, int> src2Idx;
  map<int, int> idx2Src;
  map<int, int> dest2Idx;
  map<int, int> idx2Dest;

  for (
      map<pair<int, int>, long>::const_iterator demandIt = demandMatrix.begin();
      demandIt != demandMatrix.end(); demandIt++) {
    src2Idx[demandIt->first.first] = 1;
    dest2Idx[demandIt->first.second] = 1;
  }

  map<int, int>::iterator mapIt;
  //conversion
  int numSrc = (int) src2Idx.size();
  int numDst = (int) dest2Idx.size();
  int numEdge = (int) demandMatrix.size();

  //sanity check to make sure graph
  // does not exceed capacity of hopcroft matching
  if (numSrc >= MAXN1 || numDst >= MAXN2 || numEdge > MAXM) {
    cout << "[SchedulerOptc::DemandMatchingHop] Hopcroft matching ERROR: "
         << "demand graph to large to handle! "
         << "Return without doing anything. Might result in infinite loop."
         << endl;
    return;
  }

  int nodeIdx;
  for (mapIt = src2Idx.begin(), nodeIdx = 0;
       mapIt != src2Idx.end(); mapIt++, nodeIdx++) {
    mapIt->second = nodeIdx;
    idx2Src[nodeIdx] = mapIt->first;
  }
  for (mapIt = dest2Idx.begin(), nodeIdx = 0;
       mapIt != dest2Idx.end(); mapIt++, nodeIdx++) {
    mapIt->second = nodeIdx;
    idx2Dest[nodeIdx] = mapIt->first;
  }

  int *hopc_ret_matching = (int *) malloc(numDst * sizeof(int));

  hopcroftInit(numSrc, numDst);
  for (
      map<pair<int, int>, long>::const_iterator demandIt = demandMatrix.begin();
      demandIt != demandMatrix.end(); demandIt++) {
    int srcIdx = src2Idx[demandIt->first.first];
    int destIdx = dest2Idx[demandIt->first.second];
    hopcroftAddEdge(srcIdx, destIdx);
  }
  hopcroft(hopc_ret_matching);

  for (int i = 0; i < numDst; i++) {
    if (hopc_ret_matching[i] < 0) {
      continue;
    }
    int src = idx2Src[hopc_ret_matching[i]];
    int dest = idx2Dest[i];
    mSetPtr->insert(make_pair(src, dest));
  }

  free(hopc_ret_matching);

}

// BigSliceDecomp sub-routine in Solstice (CoNext'15)
void
SchedulerOptc::BigSliceDecomp(const map<pair<int, int>, long> &demandMatrix,
                              double timeSinceCalStart) {

  if (!m_nextCircuitVector.empty() && !demandMatrix.empty()) {
    cout << "try to generate new assignments "
         << " for non-empty next circuit vector!" << endl;
  }

  cout << "now generating patterns" << endl;

  map<pair<int, int>, long> local_demand_matrix = demandMatrix;

  // Step 1: find max demnd - upperbound of threshold
  long maxDemand = -1; // used to update lower bound
  for (map<pair<int, int>, long>::const_iterator
           demandIt = local_demand_matrix.begin();
       demandIt != local_demand_matrix.end();
       demandIt++) {
    if (maxDemand < demandIt->second || maxDemand < 0) {
      maxDemand = demandIt->second;
    }
  }
  long threshold = 1;

  // max power of 2 less than max demand
  while (threshold * 2 <= maxDemand) {
    threshold = threshold * 2;
  }

  // Step 2: init matching
  set<pair<int, int>> init_matching;
  DemandMatchingHop(local_demand_matrix, &init_matching);
  unsigned long desired_matching_num = init_matching.size();

  struct timeval pattern_cal_start_time;
  struct timeval pattern_cal_end_time;

  while (local_demand_matrix.size() > 0) {

    set<pair<int, int>> *circuits = new set<pair<int, int>>();

    // STEP 1: get start time
    gettimeofday(&pattern_cal_start_time, NULL);

    // STEP 2: find a matching
    //******* threshold matching begin *********
    //                                          *
    // this couple with subtracting min_demand
    // BigSlice subroutine in Solstice
    while (circuits->size() < desired_matching_num) {
      // loop until we have a perfect matching.
      // half threshold when it is too high.
      DemandMatchingThreshold(local_demand_matrix, circuits, threshold);
      if (circuits->size() < desired_matching_num) {
        threshold /= 2;
        if (threshold < 1) {
          cout << "Error! : threshold < 1" << endl;
          exit(-1);
          // break;
          // while (circuits->size() < desired_matching_num)
        }
      }
    }
    //                                          *
    //******* threshold matching end ***********

    gettimeofday(&pattern_cal_end_time, NULL);

    // STEP 3: find circuit assignment metrics

    // STEP 3.1: demand served
    //calculate the demand served by this pattern of circuits
    //note: use max demand among the matching
    long mMinDemand = -1;

    //******* find MIN entry begin *************************
    //                                                     *
    for (set<pair<int, int>>::iterator setIt = circuits->begin();
         setIt != circuits->end(); setIt++) {
      int src = setIt->first;
      int dest = setIt->second;
      long demand = local_demand_matrix[make_pair(src, dest)];
      if (mMinDemand > demand || mMinDemand < 0) {
        mMinDemand = demand;
      }
    }
    //                                                     *
    //******* find MIN entry end ***************************

    // STEP 3.2: slot length based on demand served
    double slot_length = SecureFinishTime(mMinDemand, OPTC_BPS);

    // STEP 3.3: slot efficiency score
    double slot_score = 0.0;

    // STEP 3.4: scheduler delay
    double patternCalTime =
        secondPass(pattern_cal_end_time, pattern_cal_start_time);

    // test of pattern computation delay
    if (ZERO_COMP_TIME) {
      timeSinceCalStart = 0.0;
      patternCalTime = 0.0;
    }

    timeSinceCalStart += patternCalTime;
    double pattern_avail_time = m_currentTime + timeSinceCalStart;

    // STEP 3.5: a new assignment is generated!
    SLOT *slotp =
        new SLOT(slot_length, circuits, slot_score, pattern_avail_time);
    m_nextCircuitVector.push_back(slotp);

    if (DEBUG_LEVEL >= 10) {
      cout << __func__ << " " << circuits->size() << " matching ("
           << slot_length << ") [";
      for (set<pair<int, int>>::const_iterator
               circuit = circuits->begin();
           circuit != circuits->end();
           circuit++) {
        cout << circuit->first << "->" << circuit->second << ",";
      }
      cout << "] " << endl;
    }

    //Step 4: update demand matrix and remove zero demand
    for (set<pair<int, int>>::iterator setIt = circuits->begin();
         setIt != circuits->end(); setIt++) {
      int src = setIt->first;
      int dest = setIt->second;
      // use max demand
      local_demand_matrix[make_pair(src, dest)] -= mMinDemand;
      if (local_demand_matrix[make_pair(src, dest)] <= 0) {
        //remove entry
        // map<pair<int, int>, long>::iterator demand2RmIt;
        // demand2RmIt = demandMatrix.find(make_pair(src, dest));
        local_demand_matrix.erase(make_pair(src, dest));
      }
    }
  }
}

void
SchedulerOptc::DemandMatchingThreshold(const map<pair<int, int>,
                                                 long> &demandMatrix,
                                       set<pair<int, int>> *mSetPtr,
                                       long threshold) {

  map<pair<int, int>, long> demand_hop;
  for (
      map<pair<int, int>, long>::const_iterator demandIt = demandMatrix.begin();
      demandIt != demandMatrix.end(); demandIt++) {
    if (demandIt->second > threshold) {
      demand_hop.insert(make_pair(demandIt->first, demandIt->second));
    }
  }

  mSetPtr->clear();
  DemandMatchingHop(demand_hop, mSetPtr);

  if (DEBUG_LEVEL >= 15) {

  }
}

// takes a demand matrix as a graph
// returns the Edmond matching on the demand matrix as a circuit set
// demand should not be modified

// takes a demand matrix as a graph
// calculate Edmond (max weighted) matching on the demand matrix
// stores the matching in the set designated by the circuit set pointer
// demand matrix should not be modified
void
SchedulerOptc::DemandMatchingEdmond(const map<pair<int, int>,
                                              long> &demandMatrix,
                                    set<pair<int, int>> *mSetPtr) {
  map<int, int> src2Idx;
  map<int, int> idx2Src;
  map<int, int> dest2Idx;
  map<int, int> idx2Dest;

  map<pair<int, int>, long> edmondMatrix = demandMatrix;

  long max_trafficVol = -1;
  // scale the demand entries so as to fit into the range requirement
  for (
      map<pair<int, int>, long>::const_iterator demandIt = demandMatrix.begin();
      demandIt != demandMatrix.end(); demandIt++) {
    src2Idx[demandIt->first.first] = 1;
    dest2Idx[demandIt->first.second] = 1;
    // also find out the max traffic demand to bound
    if (max_trafficVol < demandIt->second || max_trafficVol < 0) {
      max_trafficVol = demandIt->second;
    }
  }
  map<int, int>::iterator mapIt;

  int numSrc = (int) src2Idx.size();
  int numDst = (int) dest2Idx.size();

  int nodeIdx;
  // node index starts from one
  for (mapIt = src2Idx.begin(), nodeIdx = 1;
       mapIt != src2Idx.end(); mapIt++, nodeIdx++) {
    mapIt->second = nodeIdx;
    idx2Src[nodeIdx] = mapIt->first;
  }
  // consecutive node index
  for (mapIt = dest2Idx.begin();
       mapIt != dest2Idx.end(); mapIt++, nodeIdx++) {
    mapIt->second = nodeIdx;
    idx2Dest[nodeIdx] = mapIt->first;
  }

  if (max_trafficVol > MWM_RANGE) {
    // shrink the edmondMatrix
    long double factor = (long double) MWM_RANGE / (long double) max_trafficVol;

    for (map<pair<int, int>, long>::iterator edDemandIt = edmondMatrix.begin();
         edDemandIt != edmondMatrix.end(); edDemandIt++) {
      edDemandIt->second = factor * edDemandIt->second;
    }

  } else {
    // nothing to do
  }

  Graph graph;
  int *Mate;

  graph = NewGraph(numSrc + numDst);

  for (map<pair<int, int>, long>::iterator demandIt = edmondMatrix.begin();
       demandIt != edmondMatrix.end(); demandIt++) {
    int srcIdx = src2Idx[demandIt->first.first];
    int destIdx = dest2Idx[demandIt->first.second];
    int elabel = (int) demandIt->second; // is scaled down
    AddEdge(graph, srcIdx, destIdx, elabel);
  }

  Mate = Weighted_Match(graph, 1, 1);

  //set<pair<int,int>> *circuits = new set<pair<int,int>>();

  for (int srcIdx = 1; srcIdx <= numSrc; srcIdx++) {
    // no scheduling
    if (Mate[srcIdx] == 0) continue;
    //printf("%d %d\n",srcIdx,Mate[srcIdx]);
    int src = idx2Src[srcIdx];
    int dest = idx2Dest[Mate[srcIdx]];
    mSetPtr->insert(make_pair(src, dest));
  }

  //clean up the edmond matching graph
  for (map<pair<int, int>, long>::iterator demandIt = edmondMatrix.begin();
       demandIt != edmondMatrix.end(); demandIt++) {
    int srcIdx = src2Idx[demandIt->first.first];
    int destIdx = dest2Idx[demandIt->first.second];
    Edge edge2rm = FindEdge(graph, srcIdx, destIdx);
    if (edge2rm) RemoveEdge(graph, edge2rm);
  }
  free(graph);
  free(Mate);

}

void
SchedulerOptc::SortAndSaveLargestNSlots(vector<SLOT *> &slots, int N) {
  std::sort(slots.begin(), slots.end(), slotLengthComp);
  if (N > 0 && N < (int) slots.size()) {
    vector<SLOT *> slots_to_clean(slots.begin() + N, slots.end());
    cout << "cleaning " << slots_to_clean.size() << " slots " << endl;
    deepClearSlotVector(slots_to_clean);
    slots.erase(slots.begin() + N, slots.end());
  }
}

// calculate flow rates as if the circuit is set up
// regardless of whether the real circuit is set up
void
SchedulerOptc::RateControlOptcFair(map<int, long> &optcFairRate,
                                   const long _OPTC_BPS) {

  optcFairRate.clear();

  // fair sharing in pure optical network
  map<pair<int, int>, int> circuitFlowNum;

  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {
    // for each coflow

    if ((*cfIt)->IsComplete()) {
      //such coflow has completed
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
      MapWithInc(circuitFlowNum, make_pair(src, dst), 1);
    }
  }

  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {
    // for each coflow

    if ((*cfIt)->IsComplete()) {
      //such coflow has completed
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
      int sharingFlowNum = FindWithDef(circuitFlowNum,
                                       make_pair(src, dst),
                                       1);
      long optcBps = _OPTC_BPS / sharingFlowNum;
      optcFairRate[(*fpIt)->GetFlowId()] = optcBps;
    }
  }

  if (DEBUG_LEVEL >= 5) {
    map<int, long> rates_debug = optcFairRate;
    vector<Coflow *> &coflows = m_coflowPtrVector;
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

void
SchedulerOptc::RateControlBySortedCoflows(vector<Coflow *> &coflows,
                                          map<int, long> &rates) {
  // map from flow ID to rate.
  map<int, long> expRateOptc;

  // a set of circuits that already has been taken up by other predecessor flows.
  set<pair<int, int>> circuit_has_predecessors;
  for (vector<Coflow *>::iterator cfIt = coflows.begin();
       cfIt != coflows.end(); cfIt++) {
    // for each coflow

    if ((*cfIt)->IsComplete() || 0 == (*cfIt)->GetAlpha()) {
      // such coflow has completed
      // or has zero demand
      continue;
    }
    map<pair<int, int>, long> one_coflow_allocated_demand;

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();
    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }
      std::pair<int, int> src_dst_pair((*fpIt)->GetSrc(), (*fpIt)->GetDest());
      if (ContainsKey(circuit_has_predecessors, src_dst_pair)) {
        // the circuit has been taken up by other coflows.
        continue;
      }
      MapWithInc(one_coflow_allocated_demand,
                 src_dst_pair,
                 (*fpIt)->GetBitsLeft());
    }

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      // for each flow within the coflow
      if ((*fpIt)->GetBitsLeft() <= 0) {
        //such flow has completed
        continue;
      }
      std::pair<int, int> src_dst_pair((*fpIt)->GetSrc(), (*fpIt)->GetDest());
      if (ContainsKey(circuit_has_predecessors, src_dst_pair)) {
        // the circuit has been taken up by other coflows.
        continue;
      }
      long coflow_total_demand_behind_circuit
          = MapWithDef(one_coflow_allocated_demand, src_dst_pair, (long) 0);
      if (coflow_total_demand_behind_circuit == 0) {
        cout << " interal consistency error with coflow demand " << endl;
      }
      // all flows within a coflow that shares the same circuit is going to share
      // the total bandwidth proportionally.
      long flow_rate_bitps = OPTC_BPS
          * ((double) (*fpIt)->GetBitsLeft()
              / coflow_total_demand_behind_circuit);
      MapWithInc(expRateOptc, (*fpIt)->GetFlowId(), flow_rate_bitps);
    }

    // declare circuits taken by this coflow.
    for (map<pair<int, int>, long>::iterator
             circuit_demand_iter = one_coflow_allocated_demand.begin();
         circuit_demand_iter != one_coflow_allocated_demand.end();
         circuit_demand_iter++) {
      circuit_has_predecessors.insert(circuit_demand_iter->first);
    }
  }

  rates.clear();
  rates = expRateOptc;

  if (DEBUG_LEVEL >= 5) {
    map<int, long> rates_debug = rates;
    // vector<Coflow*>& coflows = coflows;
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

void
SchedulerOptc::SchedulerAlarmPortal(double alarmTime) {

  while (!m_myTimeLine->isEmpty()) {
    Event *currentEvent = m_myTimeLine->PeekNext();
    double currentEventTime = currentEvent->GetEventTime();

    if (currentEventTime > alarmTime) {
      break;
    }

    //debug
    if (DEBUG_LEVEL >= 7) {
      cout << fixed << setw(FLOAT_TIME_WIDTH) << currentEventTime << "s "
           << "[SchedulerOptc::SchedulerAlarmPortal] "
           << " working on event type " << currentEvent->GetEventType() << endl;
    }

    bool has_flow_finished = false;
    if (m_currentTime < alarmTime) {
      has_flow_finished = Transmit(m_currentTime, alarmTime,
                                   true, // basic
                                   true, // local - should be useless
                                   FLOW_FINISH
                                       == currentEvent->GetEventType() // salvage
      );
      // validate tx constraints
      long bound =
          (alarmTime - m_currentTime) * OPTC_BPS * NUM_LINK_PER_RACK + 12;
      bool tx_meet_bound = ValidateLastTxMeetConstraints(bound);
      if (!tx_meet_bound) {
        cout << "Error: Fail to meet tx bound!!!" << endl;
        exit(-1);
      }

      m_currentTime = alarmTime;
    }

    switch (currentEvent->GetEventType()) {
      case RESCHEDULE:Schedule();
        break;
      case COFLOW_ARRIVE:CoflowArrivePreProcess();
        // clear any local traffic before calling handler.
        Transmit(m_currentTime,
                 m_currentTime, /*basic*/
                 false, /*local*/
                 true, /*salvage*/
                 false);
        // call children's coflow arrive event handler
        CoflowArriveCallBack();
        break;
      case FLOW_ARRIVE:FlowArrivePreProcess();
        // clear any local traffic before calling handler.
        Transmit(m_currentTime,
                 m_currentTime, /*basic*/
                 false, /*local*/
                 true, /*salvage*/
                 false);
        // call children's flow arrive event handler
        FlowArriveCallBack();
        break;
      case APPLY_CIRCUIT:ApplyCircuit();
        break;
      case ORDER_CIRCUIT:OrderCircuit();
        break;
      case APPLY_NEW_SCHEDULE:ApplyNewSchedule();
        break;
      case SCHEDULE_END:ScheduleEnd();
        break;
      case FLOW_FINISH:break;
      default:break;
    }
    Event *e2p = m_myTimeLine->PopNext();
    if (e2p != currentEvent) {
      cout << "[SchedulerOptc::SchedulerAlarmPortal] error: "
           << " e2p != currentEvent " << endl;
    }
    delete currentEvent;
  }

  UpdateAlarm();
}

double
SchedulerOptc::GetSlotEfficiency(set<pair<int, int>> *circuits,
                                 map<pair<int, int>, long> &demand) {
  double ret = 0.0;
  if (SORT_SLOTS_REAL_DEMAND) {
    for (set<pair<int, int>>::iterator cIt = circuits->begin();
         cIt != circuits->end(); cIt++) {
      if (demand.find(*cIt) != demand.end()) {
        ret += 1.0;
      }
    }
  }
  return ret;
}

void
SchedulerOptc::ApplyNewSchedule() {
  if (m_nextCircuitVector.empty()) {
    cout
        << "[SchedulerOptc::ApplyNewSchedule] error: No new schedule is available"
        << endl;
    return;
  }
  if (!m_circuitVector.empty()) {
    // cout << "[SchedulerOptc::ApplyNewSchedule] warning: old schedule is not ended" << endl;
    // return;
  }

  Print();

  // copy next circuit schedule to current circuit schedule
  deepClearSlotVector(m_circuitVector);
  // mem leak ?
  m_circuitVector = m_nextCircuitVector;
  m_nextCircuitVector.clear();

  // set flow rate based on m_nextElecRate & m_nextOptcRate.
  // flow finish event will be updated in OrderCircuit();
  SetFlowRate();

  m_myTimeLine->RemoveSingularEvent(SCHEDULE_END);
  m_myTimeLine->RemoveSingularEvent(ORDER_CIRCUIT);
  m_myTimeLine->RemoveSingularEvent(APPLY_CIRCUIT);

  OrderCircuit();
}

void
SchedulerOptc::deepClearSlotVector(vector<SLOT *> &v2clear) {
  for (vector<SLOT *>::iterator it = v2clear.begin();
       it != v2clear.end(); it++) {
    if (*it) {
      delete (*it)->circuits;
      delete (*it);
    }
  }
  v2clear.clear();
}

void
SchedulerOptc::ScheduleEnd() {

  if (m_circuitVector.size() > 0) {
    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerOptc::NotifyScheduleEnd] Error: "
         << " we still have available schedules. Not end!" << endl;
    return;
  }
  if (m_nextCircuitVector.size() > 0) {
    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerOptc::NotifyScheduleEnd] warmning: "
         << " new schedules are on the way. No need to reschedule. " << endl;
    return;
  }

  double checkFlows = CalcTime2FirstFlowEnd();
  if (checkFlows > 0) {
    // we have unfinished flows
    // reschedule immediately!
    Scheduler::UpdateRescheduleEvent(m_currentTime);
  }
}

bool
SchedulerOptc::OrderCircuit(void) {
  if (m_circuitVector.empty()) {
    cout << "[SchedulerOptc::OrderCircuit]"
         << " both m_circuitVector and m_nextCircuitVector are empty => no scheduling"
         << endl;

    return false;
  }

  //update circuits
  // nextActiveCircuit = circuit_to_order && current_circuits
  // circuits_switching_to = circuit_to_order - current_circuits;
  set<pair<int, int>> nextActiveCircuit;
  set<pair<int, int>> circuits_switching_to;
  set<pair<int, int>> *circuit2apply = ((*m_circuitVector.begin())->circuits);
  for (set<pair<int, int>>::iterator circuit2applyIt = circuit2apply->begin();
       circuit2applyIt != circuit2apply->end(); circuit2applyIt++) {
    if (m_activeCircuits.find(*circuit2applyIt) != m_activeCircuits.end()) {
      // the current active is also valid in the next slot
      nextActiveCircuit.insert(*circuit2applyIt);
    } else {
      circuits_switching_to.insert(*circuit2applyIt);
    }
  }
  // update active circuits set.
  m_activeCircuits.clear();
  m_activeCircuits = nextActiveCircuit;
  SetFlowPath();
  UpdateFlowFinishEvent(m_currentTime);
  // update switching circuits set.
  if (!IGNORE_WARNING_SWITCHING_CIRCUIT_SET
      && !m_circuits_switching_to.empty()) {
    cout << "Warning: m_circuits_switching_to not empty in order circuit "
         << endl;
    cout << "         Trying to order circuits in the middle of another order?"
         << endl;
    cout << "         Or applyingCircuit event removed due to new schedule?"
         << endl;
    cout
        << "          Should be Okay if you find this message around reschedule"
        << " event." << endl;
  }
  m_circuits_switching_to.clear();
  m_circuits_switching_to = circuits_switching_to;

  Event *applyCircuitEventPtr = new Event(APPLY_CIRCUIT,
                                          m_currentTime + CIRCUIT_RECONF_DELAY);
  m_myTimeLine->AddEvent(applyCircuitEventPtr);

  if (DEBUG_LEVEL >= 9) {
    cout << __func__ << " active circuits [";
    for (set<pair<int, int>>::const_iterator
             circuit_iter = m_activeCircuits.begin();
         circuit_iter != m_activeCircuits.end();
         circuit_iter++) {
      cout << circuit_iter->first << "->" << circuit_iter->second << ", ";
    }
    cout << "]" << endl;
  }

  return true; // circuits changed
}

bool
SchedulerOptc::ApplyCircuit(void) {
  if (m_circuitVector.empty()) {
    cout << fixed << setw(FLOAT_TIME_WIDTH)
         << m_currentTime << "s "
         << "[SchedulerOptc::ApplyCircuit]"
         << " Error: can't order circuit if m_circuitVector is empty!" << endl;

    return false;
  }
  m_activeCircuits.clear();

  SLOT *nextSlot = (m_circuitVector[0]);
  double slot_length = nextSlot->timeLength;
  m_activeCircuits = *(nextSlot->circuits);

  // ValidateCircuitAssignment(m_activeCircuits);

  SetFlowPath();
  UpdateFlowFinishEvent(m_currentTime);

  // update switching circuits set.
  m_circuits_switching_to.clear();

  // add order circuit event for next slot, if any
  if (m_circuitVector.size() > 1) {
    ////////////////////////////////////////////
    // use consumer/producer pattern buffer model
    double next_pattern_avail_time = m_circuitVector[1]->availTime;
    double nextOrderTime =
        max(m_currentTime + slot_length, next_pattern_avail_time);
    Event *orderCircuitEventPtr = new Event(ORDER_CIRCUIT, nextOrderTime);
    m_myTimeLine->AddEvent(orderCircuitEventPtr);
  } else if (m_circuitVector.size() <= 1) {
    // this slot is the last slot in the schedule
    // mark as end
    m_myTimeLine->RemoveSingularEvent(SCHEDULE_END);
    double scheduleEndTime = m_currentTime + nextSlot->timeLength;
    Event *scheduleEndEventPtr = new Event(SCHEDULE_END, scheduleEndTime);
    m_myTimeLine->AddEvent(scheduleEndEventPtr);

  }

  // pop out the first slot if it has been applied
  delete m_circuitVector[0]->circuits;
  delete m_circuitVector[0];
  m_circuitVector.erase(m_circuitVector.begin());

  if (DEBUG_LEVEL >= 9) {
    cout << __func__ << " active circuits [";
    for (set<pair<int, int>>::const_iterator
             circuit_iter = m_activeCircuits.begin();
         circuit_iter != m_activeCircuits.end();
         circuit_iter++) {
      cout << circuit_iter->first << "->" << circuit_iter->second << ", ";
    }
    cout << "]" << endl;
  }

  return true; // circuits changed
}

void
SchedulerOptc::SetFlowPath() {
  for (vector<Coflow *>::iterator cfIt = m_coflowPtrVector.begin();
       cfIt != m_coflowPtrVector.end(); cfIt++) {

    vector<Flow *> *flowVecPtr = (*cfIt)->GetFlows();

    for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
         fpIt != flowVecPtr->end(); fpIt++) {
      int src = (*fpIt)->GetSrc();
      int dst = (*fpIt)->GetDest();

      bool use_optcs = false;
      // only keep checking if use_optcs is false.
      for (int i = 0; i < NUM_LINK_PER_RACK && !use_optcs; i++) {
        int src_circuit = src + i * NUM_RACK;
        for (int j = 0; j < NUM_LINK_PER_RACK && !use_optcs; j++) {
          int dst_circuit = dst + j * NUM_RACK;
          if (ContainsKey(m_activeCircuits, make_pair(src_circuit,
                                                      dst_circuit))) {
            use_optcs = true;
          }
        }
      }

      if (use_optcs) {
        (*fpIt)->SetThruOptic(true);
      } else {
        (*fpIt)->SetThruOptic(false);
      }
    }
  }
}

string
SchedulerOptc::CircuitsToString(const set<pair<int, int>> &circuits) {
  if (DEBUG_LEVEL < 30) {
    return "";
  }

  std::stringstream ss;
  ss << "[";
  for (set<pair<int, int>>::const_iterator
           circuit = circuits.begin();
       circuit != circuits.end();
       circuit++) {

    int src = circuit->first;
    int dst = circuit->second;

    ss << src << "->" << dst << " , ";
  }
  ss << "]";
  return ss.str();
}

void
SchedulerOptc::CoflowArrivePreProcess() {
  EventCoflowArrive
      *coflowsArriveEvent = (EventCoflowArrive *) m_myTimeLine->PeekNext();
  if (coflowsArriveEvent->GetEventType() != COFLOW_ARRIVE) {
    cout << "[SchedulerOptc::CoflowArrivePreProcess] error: "
         << " the event type is not COFLOW_ARRIVE!" << endl;
    return;
  }

  AddCoflows(coflowsArriveEvent->m_cfpVp);
}

void
SchedulerOptc::FlowArrivePreProcess() {
  EventFlowArrive
      *flowArriveEvent = (EventFlowArrive *) m_myTimeLine->PeekNext();
  if (flowArriveEvent->GetEventType() != FLOW_ARRIVE) {
    cout << "[SchedulerOptc::FlowArrivePreProcess] error: "
         << " the event type is not FLOW_ARRIVE!" << endl;
    return;
  }

  AddFlows();
}

void
SchedulerOptc::AddCoflows(vector<Coflow *> *cfsPtr) {
  if (!cfsPtr) return;
  // check all coflows are legal
  for (vector<Coflow *>::iterator cfpIt = cfsPtr->begin();
       cfpIt != cfsPtr->end(); cfpIt++) {
    if (!(*cfpIt)->isRawCoflow()) {
      cout << "[SchedulerOptc::AddCoflows]"
           << "error: is not a raw coflow" << endl;
      continue;
    }
    m_coflowPtrVector.push_back(*cfpIt);
  }
}

void
SchedulerOptc::AddFlows() {

}