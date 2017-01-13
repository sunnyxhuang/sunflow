//
//  trafficGen.h
//  Ximulator
//
//  Created by Xin Sunny Huang on 9/21/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//

#ifndef TRAFFICGEN_H
#define TRAFFICGEN_H

#include <vector>
#include <map>
#include <set>
#include <fstream>
#include "global.h"

using namespace std;

class Flow;
class Coflow;
class JobDesc;
class Simulator;
class TrafficGenTimeLine;

class TrafficGen {
 public:
  TrafficGen();
  virtual ~TrafficGen();

  void InstallSimulator(Simulator *simPtr) { m_simPtr = simPtr; }
  /* called by simulator */
  virtual void NotifySimStart() = 0;
  virtual void TrafficGenAlarmPortal(double time) = 0;

  /* called by simulator */
  virtual void NotifyTrafficFinish(double alarmTime,
                                   vector<Coflow *> *cfpVp,
                                   vector<Flow *> *fpVp) = 0;
  virtual void NotifySimEnd() = 0;

 protected:
  Simulator *m_simPtr;
  TrafficGenTimeLine *m_myTimeLine;
  double m_currentTime;

  // internal auditing
  double m_totalCCT;
  double m_totalFCT;
  int m_totalCoflowNum;
  int m_total_accepted_coflow_num;
  int m_total_met_deadline_num;
  ifstream m_jobTraceFile;
  ofstream m_jobAuditFile;
  bool m_audit_file_title_line;

  void UpdateAlarm(); // update alarm on simulator

  friend class Simulator;
};

////////////////////////////////////////////////////
///////////// Code for FB Trace Replay   ///////////
////////////////////////////////////////////////////
class TGTraceFB : public TrafficGen {
 public:
  TGTraceFB();
  virtual ~TGTraceFB();

  /* called by simulator */
  void NotifySimStart();
  void TrafficGenAlarmPortal(double time);

  /* called by simulator */
  virtual void NotifyTrafficFinish(double alarmTime,
                                   vector<Coflow *> *cfpVp,
                                   vector<Flow *> *fpVp);
  void NotifySimEnd();

 private:
  vector<JobDesc *> m_finishJob;
  vector<JobDesc *> m_runningJob;
  vector<JobDesc *> m_readyJob;

  virtual vector<JobDesc *> ReadJobs();
  virtual void DoSubmitJob();

  void CreateFlowsWithExactSize(double time, int numMap, int numRed,
                                const vector<int> &mapLoc,
                                const vector<int> &redLoc,
                                const vector<long> &redInput,
                                vector<Flow *> &fpV);
  void CreateFlowsWithEqualSizeToSameReducer(double time,
                                             int numMap,
                                             int numRed,
                                             const vector<int> &mapLoc,
                                             const vector<int> &redLoc,
                                             const vector<long> &redInput,
                                             vector<Flow *> &fpV);
  // flow sizes will be +/- 1MB * perturb_perc%
  // perturb_perc as a percentage : perturb_perc = 10, flow sizes will be +/- 0.1 MB
  // only allow flow >= 1MB.
  void InitSeedForCoflows(int seed_for_seed);
  int GetSeedForCoflow(int coflow_id);
  void CreateFlowsWithSizePerturb(double time, int numMap, int numRed,
                                  const vector<int> &mapLoc,
                                  const vector<int> &redLoc,
                                  const vector<long> &redInput,
                                  vector<Flow *> &fpV,
                                  int perturb_perc,
                                  unsigned int rand_seed);

 protected:
  Coflow *CreateCoflowPtrFromString(double, int, int, int, string, bool, bool);
  void AnalyzeOneCoflow(Coflow *cfp, int jobid, int map, int red,
                        string &title_string, string &info_string);
  void ScheduleToAddJobs(vector<JobDesc *> &jobs);
  void KickStartReadyJobsAndNotifyScheduler();
  map<Coflow *, JobDesc *> m_coflow2job;
  vector<unsigned int> m_seed_for_coflow;
};

class JobDesc {
 public:
  JobDesc(int iid,
          double offTime,
          int numSrc,
          int numDst,
          int numFlow,
          Coflow *cfp) : m_id(iid),
                         m_offArrivalTime(offTime),
                         m_numSrc(numSrc),
                         m_numDst(numDst),
                         m_numFlow(numFlow),
                         m_coflow(cfp) {}
  ~JobDesc() {}

  int m_id; // task id
  double m_offArrivalTime; // offset arrival time
  int m_numSrc;
  int m_numDst;
  int m_numFlow;
  Coflow *m_coflow; // one coflow per task
};

///////////////////////////////////////////////////////////////////
///////////// Code for back-to-back coflow              ///////////
///////////// Replay FB trace but ignore real arrival time  ///////
///////////// so that one coflow is served at a time.     /////////
///////////////////////////////////////////////////////////////////
class TGFBOnebyOne : public TGTraceFB {
 public:
  TGFBOnebyOne();
  virtual ~TGFBOnebyOne();
 private:
  // only submit the job.
  // not reading next job.
  void DoSubmitJob();
  // add a coflow when one finishes.
  virtual void NotifyTrafficFinish(double alarmTime,
                                   vector<Coflow *> *cfpVp,
                                   vector<Flow *> *fpVp);

  // return a coflow/job at a time until the end of the trace.
  vector<JobDesc *> ReadJobs();

  void LogComptimeStat(vector<Coflow *> *cfpVp);
  ofstream m_comptimeAuditFile;
  bool m_comptime_audit_title_line;

 protected:
  // used by child class TGWindowOnebyOne
  double m_last_coflow_finish_time;
};

///////////////////////////////////////////////////////////////////////////////
///////////// Code for quick analyze of coflows ONLY                  /////////
///////////////////////////////////////////////////////////////////////////////
class TGFBAnalyzeOnly : public TGTraceFB {
 public:
  TGFBAnalyzeOnly() : TGTraceFB() {
    AnalyzeCoflow();
    m_analyze_title_line = false;
  }
  ~TGFBAnalyzeOnly() {}
 private:
  void AnalyzeCoflow();
  bool m_analyze_title_line;
};

////////////////////////////////////////////////////////////////////////////////
///////////// Code for quick analyze of network utilizaiton ONLY     ///////////
////////////////////////////////////////////////////////////////////////////////
class TGFBUtilizationOnly : public TGTraceFB {
 public:
  TGFBUtilizationOnly() : TGTraceFB() {
    m_utilization_title_line = false;
    m_analyze_title_line = false;
    AnalyzeUtilization(m_coflows);
  }
  ~TGFBUtilizationOnly() {}
 private:
  vector<Coflow *> m_coflows;
  void LoadAllCoflows(vector<Coflow *> &load_to_me);
  void AnalyzeUtilization(vector<Coflow *> &coflows);
  void CleanUpCoflows(vector<Coflow *> &clean_me);

  typedef enum UtilizationEventType_t {
    UEVENT_ARRIVE,
    UEVENT_FORGET,
    UEVENT_SAMPLE,
  } UtilizationEventType;

  struct UtilizationEvent {
    UtilizationEvent(double time, int jid, UtilizationEventType type) :
        event_time(time), job_id(jid), event_type(type) {}
    double event_time;
    int job_id;
    UtilizationEventType event_type;
  };

  struct UtilizationEventCompare {
   public:
    bool operator()(const UtilizationEvent lhs,
                    const UtilizationEvent rhs) const {
      return lhs.event_time > rhs.event_time;
    }
  };

  bool m_utilization_title_line;
  bool m_analyze_title_line;
};

#endif /*TRAFFICGEN_H*/
