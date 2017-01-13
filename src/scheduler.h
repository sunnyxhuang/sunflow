//
//  scheduler.h
//  Ximulator
//
//  Created by Xin Sunny Huang on 9/21/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//

#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <vector>
#include <fstream>
#include <queue>
#include <set>
#include <map>

using namespace std;

class Flow;
class Coflow;
class CompTimeBreakdown;
class Simulator;
class SchedulerTimeLine;

class Scheduler {
 public:
  Scheduler();
  virtual ~Scheduler();
  void InstallSimulator(Simulator *simPtr) { m_simPtr = simPtr; }
  virtual void SchedulerAlarmPortal(double) = 0;
  void NotifySimEnd();

  virtual void NotifyAddCoflows(double, vector<Coflow *> *);
  virtual void NotifyAddFlows(double);
 protected:
  Simulator *m_simPtr;
  double m_currentTime;
  SchedulerTimeLine *m_myTimeLine;
  vector<Coflow *> m_coflowPtrVector;

  map<int, long> m_nextElecRate;          /* maps flow ID to next elec rate */
  map<int, long> m_nextOptcRate;          /* maps flow ID to next optc rate */

  map<int, long> m_sBpsFree_before_workConserv;
  map<int, long> m_rBpsFree_before_workConserv;

  // map from src/dst id -> tx in bits
  map<int, long> m_validate_last_tx_src_bits;
  map<int, long> m_validate_last_tx_dst_bits;

  // for logging "utilization" based on port activity.
  ofstream m_portAuditFile;
  bool m_port_audit_title_line_last_gone;
  bool m_port_audit_title_line_bottleneck;
  bool m_port_audit_title_line_all;

  // returns true if has flow finish during transfer.
  virtual bool Transmit(double startTime,
                        double endTime,
                        bool basic,
                        bool local,
                        bool salvage);
  virtual void ScheduleToNotifyTrafficFinish(double end_time,
                                             vector<Coflow *> &coflows_done,
                                             vector<Flow *> &flows_done);
  virtual void CoflowFinishCallBack(double finishtime) = 0;
  virtual void FlowFinishCallBack(double finishTime) = 0;

  // sort coflows based on different policiess.
  void CalAlphaAndSortCoflowsInPlace(vector<Coflow *> &coflows);
  // for mosharaf's verification
  void sortCoflowOnStaticAlphaInPlace(vector<Coflow *> &coflows);
  // sort coflows based on arrival time.
  void SortCoflowsInPlaceArrival(vector<Coflow *> &coflows);

  bool ValidateLastTxMeetConstraints(long port_bound_bits);

  void SetFlowRate();

  void UpdateAlarm();
  void UpdateRescheduleEvent(double reScheduleTime);
  void UpdateFlowFinishEvent(double baseTime);

  double SecureFinishTime(long bits, long rate);
  double CalcTime2FirstFlowEnd();
  void Print(void);
 private:
  ofstream m_txAuditFile;
  // mark down circuit activities in schedulerOptc.
  // otherwise no-op.
  virtual void CircuitAuditIfNeeded(double beginTime,
                                    double endTime) {};
  virtual void WriteCircuitAuditIfNeeded(double endTime,
                                         vector<Coflow *> &coflows_done) {};
};

class SchedulerVarys : public Scheduler {
 public:
  SchedulerVarys();
  virtual ~SchedulerVarys();
  void SchedulerAlarmPortal(double currentTime);
 protected:
  // override by Varys-Deadline.
  virtual void CoflowArrive();

 private:
  virtual void Schedule(void) = 0;
  // override by Aalo.
  virtual void AddCoflows(vector<Coflow *> *cfVecPtr);

  void ApplyNewSchedule(void);
  void AddFlows();
  void FlowArrive();
  void FlowFinishCallBack(double finishTime);
  void CoflowFinishCallBack(double finishTime);

};

class SchedulerAaloImpl : public SchedulerVarys {
 public:
  SchedulerAaloImpl();
  virtual ~SchedulerAaloImpl() {}
 private:
  virtual void Schedule(void);

  virtual void AddCoflows(vector<Coflow *> *cfVecPtr);

  void RateControlAaloImpl(vector<Coflow *> &coflows,
                           vector<vector<int>> &coflow_id_queues,
                           map<int, long> &rates,
                           long LINK_RATE);
  void UpdateCoflowQueue(vector<Coflow *> &coflows,
                         vector<vector<int>> &last_coflow_id_queues,
                         map<int, Coflow *> &coflow_id_ptr_map);
  // used to maintain the stability of coflow queues.
  vector<vector<int>> m_coflow_jid_queues;
};

class SchedulerVarysPaper : public SchedulerVarys {
 public:
  SchedulerVarysPaper() : SchedulerVarys() {}
  virtual ~SchedulerVarysPaper() {}
 private:
  virtual void Schedule(void);

  // Varys described in the paper!
  // rate control based on MADD herustic (Minimum-Allocation-for-Desired-Duration)
  // and the input port based work conservation described in the paper.
  // store flow rates in rates.
  void RateControlVarysPaper(vector<Coflow *> &coflows,
                             map<int, long> &rates,
                             long LINK_RATE);

  // work conservation described in the Varys PAPER!!!
  void RateControlWorkConservationPaper(vector<Coflow *> &coflows,
                                        map<int, long> &rates,
                                        map<int, long> &sBpsFree,
                                        map<int, long> &rBpsFree,
                                        long LINK_RATE_BPS);
};

class SchedulerVarysImpl : public SchedulerVarys {
 public:
  SchedulerVarysImpl() : SchedulerVarys() {}
  virtual ~SchedulerVarysImpl() {}
 private:
  virtual void Schedule(void);

  // Varys implemented in Github!
  // rate control based on selfish coflow
  // and the work conservation for each coflow's flow.
  // store flow rates in rates.
  void RateControlVarysImpl(vector<Coflow *> &coflows,
                            map<int, long> &rates,
                            long LINK_RATE);

 protected:
  // as used by the deadline-mode varysImpl scheduler.
  // routine used RateControlVarysImpl.
  // work conservation in the Github implementation.
  void RateControlWorkConservationImpl(vector<Coflow *> &coflows,
                                       map<int, long> &rates,
                                       map<int, long> &sBpsFree,
                                       map<int, long> &rBpsFree,
                                       long LINK_RATE_BPS);
};

// The varys simulation implementation as seen in 2016 June.
class SchedulerVarysImpl201606 : public SchedulerVarys {
 public:
  SchedulerVarysImpl201606() : SchedulerVarys() {}
  virtual ~SchedulerVarysImpl201606() {}
 private:
  virtual void Schedule(void);

  // Varys implemented in Github!
  // rate control based on selfish coflow
  // and the work conservation for each coflow's flow.
  // store flow rates in rates.
  void RateControlVarysImpl201606(vector<Coflow *> &coflows,
                                  map<int, long> &rates,
                                  long LINK_RATE);

  void CheckLinkRateConstraints(const vector<Coflow *> &coflows,
                                const map<int, long> &rates,
                                long LINK_RATE_BPS);
};

class SLOT {
 public:
  SLOT(double t,
       set<pair<int, int>> *c,
       double s,
       double at) : timeLength(t),
                    circuits(c),
                    slotScore(s),
                    availTime(at) {}
  ~SLOT() {}
  double timeLength;
  set<pair<int, int>> *circuits;
  double slotScore;
  double availTime;
};

bool slotLengthComp(SLOT *s1, SLOT *s2);

typedef enum Port_Status {
  // used for src/dst port activities.
      kPortUnknown, // default
  kPortTX,      // up and tx
  kPortIdle,    // up but not tx
  kPortSolidSwitch,   // down for switching and switching for demand
  kPortFakeSwitch,    // down for switching but switching for no demand
  kPortWait,
  // used for flow activities.
      kFlowTx,
  kFlowDown,
} PortStatus;
ostream &operator<<(std::ostream &out, const PortStatus &value);

class SchedulerOptc : public Scheduler {
 public:
  SchedulerOptc();
  virtual ~SchedulerOptc();

  void SchedulerAlarmPortal(double currentTime);

 protected:

  vector<SLOT *> m_circuitVector;           /* circuit schedules */
  vector<SLOT *> m_nextCircuitVector;       /* next circuit schedules */
  // used by schedulerSunflow.
  set<pair<int, int>> m_activeCircuits;    /* current active circuits */

  double GetSlotEfficiency(set<pair<int, int>> *, map<pair<int, int>, long> &);

  void deepClearSlotVector(vector<SLOT *> &);

  void AddCoflows(vector<Coflow *> *cfVecPtr);
  void AddFlows();

  void GetDemandAll(map<pair<int, int>, long> &);

  void DemandToDoublySM(map<pair<int, int>, long> &);

  // Compared with SolsticeQuickStuff, may reduce matrix sparsity.
  void FillSquareMatrixToDSM(map<pair<int, int>, long> &,
                             map<int, long> &, map<int, long> &,
                             long &);
  // CoNExt'15 quick stuff pre-processing
  // Better way to stuff matrix
  void SolsticeQuickStuff(map<pair<int, int>, long> &demandMatrix,
                          map<int, long> &sBitFree,
                          map<int, long> &rBitFree,
                          long &bitBudget);

  void CoflowsToLayeredDemands(vector<Coflow *> &coflows,
                               vector<map<pair<int, int>, long>> &layered_demand);
  void CoflowsToLayeredDemands(vector<Coflow *> &coflows,
                               vector<map<pair<int, int>, long>> &layered_demand,
                               vector<int> &layered_coflow_ids);

  // for matrix decomposition
  void DemandMatchingHop(const map<pair<int, int>, long> &,
                         set<pair<int, int>> *);
  void DemandMatchingThreshold(const map<pair<int, int>, long> &,
                               set<pair<int, int>> *,
                               long);
  void DemandMatchingEdmond(const map<pair<int, int>, long> &demand,
                            set<pair<int, int>> *circuit_set);

  void HopDecomp(map<pair<int, int>, long> &, double);
  void BigSliceDecomp(const map<pair<int, int>, long> &, double);

  // post decomposition procedures.
  // sort, in place, all generated slots in slots,
  // and save only the largest N slots, dumping the rest small ones.
  void SortAndSaveLargestNSlots(vector<SLOT *> &slots, int N);

  // Rate assignment based on pure optical fair sharing
  // when a circuit is set up,
  // all flows on the circuit have fair share of link bandwidth
  void RateControlOptcFair(map<int, long> &optcFairRate,
                           const long _OPTC_BPS);

  // Rate assigned in the order of coflows.
  // Pioritized coflow non-blocked by less prioritized coflow on a circuit.
  void RateControlBySortedCoflows(vector<Coflow *> &coflows,
                                  map<int, long> &rates);

  // tool function
  void SetFlowPath();

  // unbox coflows ptr
  //  and call children's CoflowArrive/FlowArrive
  void CoflowArrivePreProcess();
  void FlowArrivePreProcess();

  virtual void ApplyNewSchedule(void);

  // return true if circuit is changed
  virtual bool OrderCircuit(void);

  // return false if circuit is NOT changed
  virtual bool ApplyCircuit(void);

  virtual void ScheduleEnd();

  //debug
  // validate circuits is a matching.
  bool ValidateCircuitAssignment(const set<pair<int, int>> &);
  // validate a demand matrix to be a square doubly stochastic matrix
  void ValidateDSM(map<pair<int, int>, long> &);

  // defined in scheduler.
  // used in scheduler::Transmit.
  // mark down circuit activities.
  // in src/dst/flow activity.
  virtual void CircuitAuditIfNeeded(double beginTime,
                                    double endTime);
  virtual PortStatus CheckPortStatus(int port,
                                     bool is_src, // should be useless
                                     const pair<int, int> &circuit_if_active,
                                     const pair<int, int> &circuit_if_switching,
                                     const map<pair<int, int>, long> &demand);
  virtual void WriteCircuitAuditIfNeeded(double endTime,
                                         vector<Coflow *> &coflows_done);

  friend ostream &operator<<(std::ostream &out, const PortStatus &value);

  class PortActivity {
   public:
    PortActivity(double start,
                 double end,
                 PortStatus status) :
        m_start_time(start),
        m_end_time(end),
        m_status(status) { m_peer = -1; };
    double m_start_time;
    double m_end_time;
    PortStatus m_status;
    // valid only when status is kPortBusy.
    // -1 by default, re-init if needed.
    int m_peer;
  };
  struct LatestEndTimeFirst {
    bool operator()(const PortActivity &lhs,
                    const PortActivity &rhs) const {
      return lhs.m_end_time < rhs.m_end_time;
    }
  };

  // last activity at the front
  typedef std::priority_queue<PortActivity,
                              std::vector<PortActivity>,
                              LatestEndTimeFirst> PortActivityQueue;
  typedef map<int, PortActivityQueue> PortActivityMap;
  PortActivityMap m_src_activity;
  PortActivityMap m_dst_activity;

  // either one is valid.
  int m_max_optimal_workspan_src;
  int m_max_optimal_workspan_dst;
  // activity of port with max optimal-workspan
  PortActivityQueue m_max_optimal_workspan_port_activity;

  // set of circuits ordered but not yet applied.
  // MUST be updated open OrderCircuit and ApplyCircuit.
  set<pair<int, int>> m_circuits_switching_to;

  string CircuitsToString(const set<pair<int, int>> &circuits);

 private:

  // implemented in different scheduler
  virtual void Schedule(void) = 0;
  virtual void CoflowArriveCallBack(void)=0;
  virtual void FlowArriveCallBack()=0;
  virtual void CoflowFinishCallBack(double finishtime)=0;
  virtual void FlowFinishCallBack(double finishtime)=0;


  ofstream m_circuitAuditFile;

  /* data */
};

class SchedulerSunflow : public SchedulerOptc {
 public:
  SchedulerSunflow();
  virtual ~SchedulerSunflow();
 private:
  virtual void Schedule(void);
  virtual void CoflowArriveCallBack(void);
  virtual void FlowArriveCallBack();
  virtual void CoflowFinishCallBack(double finishtime);
  virtual void FlowFinishCallBack(double finishTime);

  // override :
  virtual void ApplyNewSchedule(void);
  // override : return true if circuit is changed
  virtual bool OrderCircuit(void);
  // override : return false if circuit is NOT changed
  virtual bool ApplyCircuit(void);

  // call sunflow solver to compute reservations.
  virtual void ComputeReservations(vector<Coflow *> &coflows);
  // return the absolute (virtual) time of reservations.
  double GetNextReservationTimeAndCircuits(set<pair<int,
                                                    int>> &circuits_to_set);

  void ExcludeConflictingCircuits(const set<pair<int, int>> &current_circuits,
                                  const set<pair<int, int>> &circuits_to_set,
                                  set<pair<int, int>> &remaining_circuits);
  void SetMinus(const set<pair<int, int>> &lhs,
                const set<pair<int, int>> &rhs,
                set<pair<int, int>> &lhs_minus_rhs);
  // list of circuit to become valid events;
  vector<pair<double, set<pair<int, int>>>> circuits_to_apply_list;

  // resource leakage check.
  void FindLeakedCircuits(const vector<Coflow *> &coflows,
                          const set<pair<int, int>> &occupied_circuits,
                          set<pair<int, int>> &leaked_circuits);
  void FindCircuitsWithDemand(const vector<Coflow *> &coflows,
                              set<pair<int, int>> &circuits_with_demand);
  void MarkPotentialResourceLeak(double alarm_time);
  bool ValidateNoResourceLeak(const vector<Coflow *> &coflows,
                              const set<pair<int, int>> &next_active_circuits,
                              double tolerance);
  double GetLeakCheckTolerance();
  // list of (m_extra_circuits_since_time, m_extra_circuits)
  vector<pair<double, set<pair<int, int>>>> m_leaked_circuits_profile;

 protected:
  typedef std::pair<double, map<pair<int, int>, double>> Reservation;
  class ReservationComp {
   public:
    bool operator()(Reservation lhs, Reservation rhs) {
      return lhs.first > rhs.first;
    }
  };
  typedef priority_queue<Reservation, vector<Reservation>, ReservationComp>
      CircuitReservations;
  CircuitReservations m_reservations;
  // circuit to expire time.
  // by default, ALL circuits expire => useless.
  map<pair<int, int>, double> m_circuit_to_expired_time;

  double m_last_schedule_valid_for;
};

class SchedulerSolstice : public SchedulerOptc {
 public:
  SchedulerSolstice();
  virtual ~SchedulerSolstice();
 private:
  // different behaviors in case of different events
  // pure virtual functions declared in schedulerOptc base class
  virtual void Schedule(void);
  virtual void CoflowArriveCallBack(void);
  virtual void FlowArriveCallBack();
  virtual void CoflowFinishCallBack(double finishtime);
  virtual void FlowFinishCallBack(double finishTime);
};

class SchedulerBvN : public SchedulerOptc {
 public:
  SchedulerBvN();
  virtual ~SchedulerBvN();
 private:
  // different behaviors in case of different events
  // pure virtual functions declared in schedulerOptc base class
  virtual void Schedule(void);
  virtual void CoflowArriveCallBack(void);
  virtual void FlowArriveCallBack();
  virtual void CoflowFinishCallBack(double finishtime);
  virtual void FlowFinishCallBack(double finishTime);
};

class SchedulerTMS : public SchedulerOptc {
 public:
  SchedulerTMS();
  virtual ~SchedulerTMS();
 private:

  // different behaviors in case of different events

  // pure virtual functions declared in schedulerOptc base class
  virtual void Schedule(void);
  virtual void CoflowArriveCallBack(void);
  virtual void FlowArriveCallBack();
  virtual void CoflowFinishCallBack(double finishtime);
  virtual void FlowFinishCallBack(double finishTime);

 protected:
  // takes in a demand matrix
  // convert to a doubly stochastic matrix that
  // each col/row sums to resulting_sum
  void Sinkhorn(map<pair<int, int>, long> &, long resulting_sum);
};

class SchedulerEdmond : public SchedulerOptc {
 public:
  SchedulerEdmond();
  virtual ~SchedulerEdmond();
 private:

  // different behaviors in case of different events

  // pure virtual functions declared in schedulerOptc base class
  virtual void Schedule(void);
  virtual void CoflowArriveCallBack(void);
  virtual void FlowArriveCallBack();
  virtual void CoflowFinishCallBack(double finishtime);
  virtual void FlowFinishCallBack(double finishTime);

  // perform edmond matching on demand matrix, return the most
  // demanding circuit set.
  void GenerateOneEdmondAssignment(map<pair<int, int>, long> &demandMatrix,
                                   double timeSinceCalStart);
};

#endif /*SCHEDULER_H*/
