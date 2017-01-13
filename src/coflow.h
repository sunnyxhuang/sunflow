//
//  coflow.h
//  Ximulator
//
//  Created by Xin Sunny Huang on 9/21/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//

#ifndef COFLOW_H
#define COFLOW_H

#include <vector>
#include "util.h"

#include <map>
#include <iostream>

using namespace std;

class Flow {
 public:
  Flow(double startTime, int src, int dest, long sizeInByte);
  ~Flow();
  int GetFlowId() { return m_flowId; }
  double GetStartTime() { return m_startTime; }
  double GetEndTime() { return m_endTime; }
  void SetEndTime(double endtime) { m_endTime = endtime; }
  int GetSrc() { return m_src; }
  int GetDest() { return m_dest; }
  long GetSizeInBit() { return m_sizeInBit; }
  long GetBitsLeft() { return m_bitsLeft; }
  long GetElecRate(void) { return m_elecBps; }
  long GetOptcRate(void) { return m_optcBps; }
  void SetRate(long elecBps, long optcBps);
  void SetThruOptic(bool thruOptc) { m_thruOptic = thruOptc; }
  bool isThruOptic() { return m_thruOptic; }
  bool isRawFlow();
  long Transmit(double startTime, double endTime);// return bitsLeft
  long TxSalvage();// clear flow if smaller than threshold.
  long TxLocal();  // clear flow if local rack; return bitsLeft otherwise
  string toString();
 private:
  static int s_flowIdTracker;
  int m_flowId;
  int m_src;
  int m_dest;
  bool m_thruOptic;
  long m_sizeInBit;
  long m_bitsLeft;
  double m_startTime;
  double m_endTime;
  long m_elecBps;
  long m_optcBps;
};

class CompTimeBreakdown {
 public:
  CompTimeBreakdown() {
    m_time_shuffle_random = 0.0;
    m_time_sort = 0.0;
    m_time_reserve = 0.0;
    m_time_stuff = 0.0;
    m_time_slice = 0.0;
  }
  double GetSunflowTotalTime() {
    return m_time_shuffle_random
        + m_time_sort
        + m_time_reserve;
  }
  double GetSolsticeTotalTime() {
    return m_time_stuff + m_time_slice;
  }
  double GetVectorAvgTime() {
    if (m_time_vector.empty()) return 0.0;
    double sum = 0.0;
    for (vector<double>::const_iterator
             t = m_time_vector.begin();
         t != m_time_vector.end();
         t++) {
      sum += *t;
    }
    return sum / m_time_vector.size(); // avg
  }
  double GetVectorMaxTime() {
    double max = -1.0;
    for (vector<double>::const_iterator
             t = m_time_vector.begin();
         t != m_time_vector.end();
         t++) {
      if (max < *t || max < 0) {
        max = *t;
      }
    }
    return max;
  }
  double GetVectorMinTime() {
    double min = -1.0;
    for (vector<double>::const_iterator
             t = m_time_vector.begin();
         t != m_time_vector.end();
         t++) {
      if (min > *t || min < 0) {
        min = *t;
      }
    }
    return min;
  }
  // for sunflow
  double m_time_shuffle_random;
  double m_time_sort;
  double m_time_reserve;
  // for solstice
  double m_time_stuff;
  double m_time_slice;
  // for other decomposition-based schedulers.
  vector<double> m_time_vector;
};

class Coflow {
 public:
  Coflow(double startTime, int totalFlows);
  ~Coflow();
  int GetCoflowId() { return m_coflowId; }
  int GetJobId() { return m_job_id; }
  void SetJobId(int job_id) { m_job_id = job_id; }
  double GetEndTime() { return m_endTime; }
  void SetEndTime(double endTime) { m_endTime = endTime; }
  void AddFlow(Flow *f);

  long GetMaxPortLoadInBits();
  double GetMaxOptimalWorkSpanInSeconds();
  long GetLoadOnMaxOptimalWorkSpanInBits();
  void GetPortOnMaxOptimalWorkSpan(int &src, int &dst);
  long GetLoadOnPortInBits(int src, int dst);
  double GetOptimalWorkSpanOnPortInSeconds(int src, int dst);
  double GetSizeInByte() { return m_coflow_size_in_bytes; }

  void AddTxBit(long bit_sent) {
    m_coflow_sent_bytes += (double) bit_sent / 8.0;
  };
  double GetSentByte() { return m_coflow_sent_bytes; };

  void SetDeadlineSec(double d) { m_deadline_duration = d; }
  double GetDeadlineSec() { return m_deadline_duration; }

  bool IsRejected() { return m_is_rejected; }
  void SetRejected() { m_is_rejected = true; }

  void SetCompTime(CompTimeBreakdown comptime) { m_comp_time = comptime; }
  void AddTimeToVector(double one_slot_computation_time) {
    m_comp_time.m_time_vector.push_back(one_slot_computation_time);
  }
  CompTimeBreakdown GetCompTime() { return m_comp_time; }

  void SetStaticAlpha(long alpha) { m_static_alpha = alpha; }
  long GetStaticAlpha() { return m_static_alpha; }

  long CalcAlpha();
  long GetAlpha() { return m_alpha; }

  double CalcBeta();
  double GetBeta() { return m_beta; }

  double CalcAlphaOnline(map<int, long> &sBpsFree,
                         map<int, long> &rBpsFree,
                         long LINK_RATE_BPS);
  double GetAlphaOnline() { return m_online_alpha; }

  double GetStartTime() { return m_startTime; }
  bool IsComplete() { return m_nFlowsCompleted >= m_nTotalFlows; }
  bool IsFlowsAddedComplete() { return m_nFlowsCompleted >= m_nFlows; }

  vector<Flow *> *GetFlows(void) { return m_flowVector; }

  bool NumFlowFinishInc();
  bool isRawCoflow();

  void Print();
  string toString();

 private:
  vector<Flow *> *m_flowVector;
  static int s_coflowIdTracker;

  int m_job_id;
  int m_coflowId;
  int m_nFlows;
  int m_nFlowsCompleted;
  int m_nTotalFlows;

  // Expected CCT based on packet switching.
  long m_static_alpha;
  long m_alpha;
  double m_online_alpha;

  // Expected CCT based on circuit switching, used for sunflow.
  double m_beta;

  double m_startTime; /*earliest start time of flow contained*/
  double m_endTime;
  // deadline duration after arrival time.
  double m_deadline_duration;
  bool m_is_rejected;

  // coflow profile.
  map<int, long> m_src_bits;
  map<int, long> m_dst_bits;

  map<int, double> m_src_time;
  map<int, double> m_dst_time;

  double m_coflow_size_in_bytes;

  // aalo
  double m_coflow_sent_bytes;

  CompTimeBreakdown m_comp_time;

  /* data */
};

bool coflowCompAlpha(Coflow *l, Coflow *r);
bool coflowCompBeta(Coflow *l, Coflow *r);
bool coflowCompArrival(Coflow *l, Coflow *r);
#endif /*COFLOW_H*/
