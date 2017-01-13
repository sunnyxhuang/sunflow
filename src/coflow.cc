//
//  coflow.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 9/21/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//

#include "coflow.h"
#include "global.h"

#include <ctgmath>
#include <iomanip> // setw

using namespace std;

int Flow::s_flowIdTracker = 0;

Flow::Flow(double startTime, int src, int dest, long sizeInByte) {
  m_flowId = s_flowIdTracker++;
  m_startTime = startTime;
  m_endTime = INVALID_TIME;
  m_src = src;
  m_dest = dest;

  m_thruOptic = false;
  m_sizeInBit = sizeInByte * 8 * TRAFFIC_SIZE_INFLATE;
  if (m_sizeInBit < 1000000 * 8) {
    m_sizeInBit = 1000000 * 8; // all flows >= 1MB;
  }
  m_bitsLeft = m_sizeInBit;
  m_elecBps = 0.0;
  m_optcBps = 0.0;

}

Flow::~Flow() {
  //cout << "Flow destructor called." << endl;
  /*do nothing */
}

void Flow::SetRate(long elecBps, long optcBps) {
  m_elecBps = elecBps;
  m_optcBps = optcBps;
}

bool
Flow::isRawFlow() {
  if (m_bitsLeft != m_sizeInBit
      || m_elecBps != 0
      || m_optcBps != 0) {
    return false;
  }
  return true;
}

// return bitsLeft
// Require endTime > startTime
long Flow::Transmit(double startTime, double endTime) {
  long tBits = 0;
  long bps = 0;
  if (m_thruOptic) {
    tBits = ceil((endTime - startTime) * m_optcBps);
    bps = m_optcBps;
  } else {
    tBits = ceil((endTime - startTime) * m_elecBps);
    bps = m_elecBps;
  }
  long tBitsValid = tBits < m_bitsLeft ? tBits : m_bitsLeft;
  m_bitsLeft -= tBitsValid;

  //debug msg
  if (DEBUG_LEVEL >= 35) {
    cout << string(FLOAT_TIME_WIDTH + 2, ' ')
         /* make room for time in other lines */
         << "[Flow::Transmit] Flow [" << m_flowId << "] "
         << m_src << "->" << m_dest
         << " TX " << tBitsValid << " bits in"
         << "(" << startTime << ", " << endTime << ")s "
         << "rate " << bps << " bps "
         << m_bitsLeft << " bits_left" << endl;
  }

  return m_bitsLeft;
}

long
Flow::TxLocal() {
  if (m_src == m_dest) {
    m_bitsLeft = 0;
  }
  return m_bitsLeft;
}

long
Flow::TxSalvage() {
  long bps = 0;
  if (m_thruOptic) {
    bps = m_optcBps;
  } else {
    bps = m_elecBps;
  }
  if (bps > 0 && m_bitsLeft <= 10) {
    m_bitsLeft = 0;
  }
  return m_bitsLeft;
}

string
Flow::toString() {
  std::stringstream ss;
  ss << "Flow-[" << m_flowId << ", "
     << m_src << "->" << m_dest << " "
     << m_bitsLeft << " bits] ";
  return ss.str();
}

void
Coflow::Print() {
  for (vector<Flow *>::iterator it = m_flowVector->begin();
       it != m_flowVector->end(); it++) {
    if ((*it)->GetBitsLeft() <= 0) {
      continue;
    }
    cout << "            " /* make room for time in other lines */
         << "[Coflow::Print] Flow Id " << (*it)->GetFlowId() << " "
         << (*it)->GetSrc() << "=>" << (*it)->GetDest() << " "
         << (*it)->GetBitsLeft() << " bits_left "
         << (*it)->GetElecRate() << " bps_elec "
         << (*it)->GetOptcRate() << " bps_optc " << endl;
  }
  cout << string(FLOAT_TIME_WIDTH + 2,
                 ' ') /* make room for time in other lines */
       << "[Coflow::Print] CalcAlpha = " << CalcAlpha() << endl;
}

string
Coflow::toString() {
  std::stringstream ss;
  ss << "Coflow-[" << m_job_id << ", "
     << m_startTime << "+"
     << GetMaxOptimalWorkSpanInSeconds() << ", "
     << "a = " << GetAlpha() << ", "
     << setw(5) << m_nFlowsCompleted << " / "
     << setw(5) << m_nTotalFlows << " flows] ";
  return ss.str();
}

int Coflow::s_coflowIdTracker = 0;

Coflow::Coflow(double startTime, int totalFlows) {
  m_job_id = -1;
  m_coflowId = s_coflowIdTracker++;
  m_startTime = startTime;
  m_nTotalFlows = totalFlows;
  m_nFlows = 0;
  m_nFlowsCompleted = 0;
  //m_finalFlowAdded = false;
  m_flowVector = new vector<Flow *>;
  m_alpha = -1;
  m_beta = -1.0;
  m_coflow_size_in_bytes = 0.0;
  m_coflow_sent_bytes = 0.0;
  //
  m_endTime = INVALID_TIME;
  m_deadline_duration = 0.0;
  m_is_rejected = false;
}

Coflow::~Coflow() {
//    cout << "Coflow destructor called." << endl;
  //delete all flows within the coflow
  for (vector<Flow *>::iterator it = m_flowVector->begin();
       it != m_flowVector->end(); it++) {
    if (*it) {
      delete (*it);
    }
  }
  delete m_flowVector;
}

bool
Coflow::NumFlowFinishInc() {
  if (m_nFlowsCompleted < m_nFlows) {
    m_nFlowsCompleted++;
    return true;
  }
  return false;
}

bool
Coflow::isRawCoflow() {
 
  if (!m_flowVector) {
    return false;
  }
  for (vector<Flow *>::iterator fpIt = m_flowVector->begin();
       fpIt != m_flowVector->end(); fpIt++) {
    if (!(*fpIt)->isRawFlow()) {
      return false;
    }
  }
  return true;
}

void
Coflow::AddFlow(Flow *fp) {
  if (!fp) return;
  m_flowVector->push_back(fp);
  m_nFlows++;

  // do some coflow profiling.
  // by default, use 1G to estimate.
  if (fp->GetSrc() != fp->GetDest()) {
    double time = fp->GetSizeInBit() / (double) DEFAULT_LINK_RATE_BPS
        + CIRCUIT_RECONF_DELAY;
    MapWithInc(m_src_time, fp->GetSrc(), time);
    MapWithInc(m_dst_time, fp->GetDest(), time);
    long bits = fp->GetSizeInBit();
    MapWithInc(m_src_bits, fp->GetSrc(), bits);
    MapWithInc(m_dst_bits, fp->GetDest(), bits);
    m_coflow_size_in_bytes += (bits / 8.0);
  }
}

long
Coflow::GetMaxPortLoadInBits() {
  double max_src = MaxMap(m_src_bits);
  double max_dst = MaxMap(m_dst_bits);
  return max_src > max_dst ? max_src : max_dst;
}

double
Coflow::GetMaxOptimalWorkSpanInSeconds() {
  double max_src = MaxMap(m_src_time);
  double max_dst = MaxMap(m_dst_time);
  return max_src > max_dst ? max_src : max_dst;
}

//return payload, in bits, on the src/dst port
// with the max optimal work span.
long
Coflow::GetLoadOnMaxOptimalWorkSpanInBits() {
  int src = -1;
  int dst = -1;
  GetPortOnMaxOptimalWorkSpan(src, dst);
  if ((src <= -1) != (dst <= -1)) {
    // solid pair : either one is valid.
    return GetLoadOnPortInBits(src, dst);
  }
  cout << __func__ << ": port pair invalid" << endl;
  return -1;
}

//return port with max optimal-work-span
void
Coflow::GetPortOnMaxOptimalWorkSpan(int &src, int &dst) {
  double max_src = MaxMap(m_src_time);
  double max_dst = MaxMap(m_dst_time);
  if (max_src > max_dst) {
    // max optimal-work-span on src.
    src = Key2MaxPositiveMap(m_src_time);
    dst = -1;
    //    cout << "max optimal work span on src " << src << " time " << max_src
    //         << endl;
  } else {
    // max optimal-work-span on dst.
    src = -1;
    dst = Key2MaxPositiveMap(m_dst_time);
    //    cout << "max optimal work span on dst " << dst << " time " << max_dst
    //         << endl;
  }
}

long
Coflow::GetLoadOnPortInBits(int src, int dst) {
  if ((src <= -1) != (dst <= -1)) {
    // solid pair : either one is valid.
    return dst <= -1 ?
           MapWithDef(m_src_bits, src, (long) 0) : MapWithDef(m_dst_bits,
                                                              dst,
                                                              (long) 0);;
  }
  cout << __func__ << ": port pair invalid" << endl;
  return -1;
}

double
Coflow::GetOptimalWorkSpanOnPortInSeconds(int src, int dst) {
  if ((src <= -1) != (dst <= -1)) {
    // solid pair : either one is valid.
    return dst <= -1 ?
           MapWithDef(m_src_time, src, 0.0) : MapWithDef(m_dst_time, dst, 0.0);;
  }
  cout << __func__ << ": port pair invalid" << endl;
  return -1;
}

// the higher, the less prioritized
long
Coflow::CalcAlpha() {
  long maxPortSumBits = 0;
  map<int, long> sBits;
  map<int, long> rBits;
  for (vector<Flow *>::iterator it = m_flowVector->begin();
       it != m_flowVector->end(); it++) {
    MapWithInc(sBits, (*it)->GetSrc(), (*it)->GetBitsLeft());
    MapWithInc(rBits, (*it)->GetDest(), (*it)->GetBitsLeft());
  }

  long sBitsMax = MaxMap(sBits);
  long rBitsMax = MaxMap(rBits);

  maxPortSumBits = sBitsMax > rBitsMax ? sBitsMax : rBitsMax;

  m_alpha = maxPortSumBits;
  return maxPortSumBits;
}

double
Coflow::CalcBeta() {
  m_beta = 0.0;

  map<int, double> src_time;
  map<int, double> dst_time;

  int num_active_flows = 0;
  for (vector<Flow *>::iterator
           it = m_flowVector->begin();
       it != m_flowVector->end();
       it++) {
    if ((*it)->GetBitsLeft() <= 0) {
      // skip finished flows
      continue;
    }
    // this flow is yet to be finished.
    double time = (*it)->GetBitsLeft() / (double) DEFAULT_LINK_RATE_BPS
        + CIRCUIT_RECONF_DELAY;
    MapWithInc(src_time, (*it)->GetSrc(), time);
    MapWithInc(dst_time, (*it)->GetDest(), time);
    num_active_flows++;
  }

  double max_src = MaxMap(src_time);
  double max_dst = MaxMap(dst_time);
  //  1st term : sunflow 2X optimality.
  double tx_delay = 2 * max(max_src, max_dst);
  //  2nd term :
  double preemption_delay = num_active_flows * CIRCUIT_RECONF_DELAY;

  m_beta = tx_delay + preemption_delay;
  return m_beta;
}

// calc min completion time for a coflow
// and store it as the alpha
// based on current bandwidth remained.
double
Coflow::CalcAlphaOnline(map<int, long> &sBpsFree,
                        map<int, long> &rBpsFree,
                        long LINK_RATE_BPS) {

  map<int, long> sBits;
  map<int, long> rBits;
  for (vector<Flow *>::iterator it = m_flowVector->begin();
       it != m_flowVector->end(); it++) {
    MapWithInc(sBits, (*it)->GetSrc(), (*it)->GetBitsLeft());
    MapWithInc(rBits, (*it)->GetDest(), (*it)->GetBitsLeft());
  }

  map<int, double> sAlphaSec;
  map<int, double> rAlphaSec;

  for (map<int, long>::const_iterator sIter = sBits.begin();
       sIter != sBits.end(); sIter++) {
    long sBps = MapWithDef(sBpsFree, sIter->first, LINK_RATE_BPS);
    if (sBps <= 0) {
      // debug
      // cout << " src " << sIter->first << " no remaining bandwidth" << endl;
      m_online_alpha = CF_DEAD_ALPHA_SIGN;
      return CF_DEAD_ALPHA_SIGN;
    }
    sAlphaSec[sIter->first] = sIter->second / (double) sBps;
    // cout << " sAlphaSec " << sAlphaSec[sIter->first] << endl;
  }

  for (map<int, long>::const_iterator rIter = rBits.begin();
       rIter != rBits.end(); rIter++) {
    long rBps = MapWithDef(rBpsFree, rIter->first, LINK_RATE_BPS);
    if (rBps <= 0) {
      // debug
      // cout << " dst " << rIter->first << " no remaining bandwidth" << endl;
      m_online_alpha = CF_DEAD_ALPHA_SIGN;
      return CF_DEAD_ALPHA_SIGN;
    }
    rAlphaSec[rIter->first] = rIter->second / (double) rBps;
    // cout << " rAlphaSec "  << rAlphaSec[rIter->first] << endl;
  }

  double sAlphaSecMax = MaxMap(sAlphaSec);
  double rAlphaSecMax = MaxMap(rAlphaSec);

  m_online_alpha = sAlphaSecMax > rAlphaSecMax ? sAlphaSecMax : rAlphaSecMax;

  // if too long ...
  if (m_online_alpha > ONLINE_ALPHA_CUTOFF) {
    m_online_alpha = CF_DEAD_ALPHA_SIGN;
    return CF_DEAD_ALPHA_SIGN;
  }

  return m_online_alpha;
}

bool coflowCompAlpha(Coflow *l, Coflow *r) {
  return l->GetAlpha() <= r->GetAlpha();
}

bool coflowCompBeta(Coflow *l, Coflow *r) {
  return l->GetBeta() <= r->GetBeta();
}

bool coflowCompArrival(Coflow *l, Coflow *r) {
  return l->GetStartTime() <= r->GetStartTime();
}