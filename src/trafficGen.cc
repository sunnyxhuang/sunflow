//
//  trafficGen.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 9/21/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//
#include <iomanip>
#include <algorithm>
#include <cfloat>
#include <ctgmath>
#include <numeric> //dbg
#include <queue> // for utilization analysis.
#include <random>
#include <sstream>

#include "util.h"
#include "global.h"
#include "trafficGen.h"
#include "scheduler.h"
#include "coflow.h"
#include "events.h"

using namespace std;

TrafficGen::TrafficGen() {
  m_currentTime = 0;
  m_simPtr = NULL;
  m_myTimeLine = new TrafficGenTimeLine;

  m_totalCCT = 0;
  m_totalFCT = 0;

  m_totalCoflowNum = 0;
  m_total_accepted_coflow_num = 0;
  m_total_met_deadline_num = 0;
  m_audit_file_title_line = false;
}

TrafficGen::~TrafficGen() {
  delete m_myTimeLine;

  cout << "Done" << " ";
  cout << m_total_met_deadline_num
       << "/" << m_total_accepted_coflow_num
       << "/" << m_totalCoflowNum << " ";
  cout << "" << m_totalCCT << " ";
  cout << "" << m_totalFCT << " ";
  cout << endl;

  if (m_jobTraceFile.is_open() && m_jobAuditFile.is_open()) {
    m_jobTraceFile.close();
    //trace files
    m_jobAuditFile << "Done" << " ";
    m_jobAuditFile << m_total_met_deadline_num
                   << "/" << m_total_accepted_coflow_num
                   << "/" << m_totalCoflowNum << " ";
    m_jobAuditFile << "" << m_totalCCT << " ";
    m_jobAuditFile << "" << m_totalFCT << " ";

    if (SUNFLOW_SHUFFLE_RANDOM
        || SUNFLOW_SORT_DEMAND) {
      // sunflow only.
      // no shuffle info by default.
      m_jobAuditFile << (SUNFLOW_SHUFFLE_RANDOM ? "ShuffleR" : "********");
      m_jobAuditFile << "-";
      m_jobAuditFile << (SUNFLOW_SORT_DEMAND ? "Sort" : "****");
      m_jobAuditFile << " ";
    }

    m_jobAuditFile << endl;
    m_jobAuditFile.close();
  }
}

void
TrafficGen::UpdateAlarm() {
  if (!m_myTimeLine->isEmpty()) {
    Event *nextEvent = m_myTimeLine->PeekNext();
    double nextTime = nextEvent->GetEventTime();
    Event *ep = new Event(ALARM_TRAFFIC, nextTime);
    m_simPtr->UpdateTrafficAlarm(ep);
  }
}

////////////////////////////////////////////////////
///////////// Code for FB Trace Replay   ///////////
////////////////////////////////////////////////////


TGTraceFB::TGTraceFB() {

  m_readyJob = vector<JobDesc *>();
  m_runningJob = vector<JobDesc *>();

  m_coflow2job = map<Coflow *, JobDesc *>();

  m_jobTraceFile.open(TRAFFIC_TRACE_FILE_NAME);
  if (!m_jobTraceFile.is_open()) {
    cout << "Error: unable to open file "
         << TRAFFIC_TRACE_FILE_NAME << endl;
    cout << "Now terminate the program" << endl;
    exit(-1);
  }
  // trace files
  m_jobAuditFile.open(TRAFFIC_AUDIT_FILE_NAME);
  if (!m_jobAuditFile.is_open()) {
    cout << "Error: unable to open file "
         << TRAFFIC_AUDIT_FILE_NAME << endl;
    cout << "Now terminate the program" << endl;
    exit(-1);
  }

  int seed_for_perturb_seed = 13;
  InitSeedForCoflows(seed_for_perturb_seed);
}

TGTraceFB::~TGTraceFB() {
  if (!m_runningJob.empty()) {
    //error
    cout << "[TGTraceFB::~TGTraceFB()] Error:"
         << "has unfinished traffic at the end of simulation!"
         << endl;
  }
}

/* called by simulator */
void
TGTraceFB::NotifySimStart() {
  vector<JobDesc *> jobs2add = ReadJobs();
  ScheduleToAddJobs(jobs2add);
  UpdateAlarm();
}

void
TGTraceFB::TrafficGenAlarmPortal(double alarmTime) {

  while (!m_myTimeLine->isEmpty()) {
    Event *currentEvent = m_myTimeLine->PeekNext();
    double currentEventTime = currentEvent->GetEventTime();

    if (currentEventTime > alarmTime) {
      // break if we have over run the clock
      break;
    }
    // it seems everything is ok
    m_currentTime = alarmTime;

    cout << fixed << setw(FLOAT_TIME_WIDTH) << currentEventTime << "s "
         << "[TGTraceFB::TrafficGenAlarmPortal] "
         << " working on event type " << currentEvent->GetEventType() << endl;

    //currentEvent->CallBack();
    switch (currentEvent->GetEventType()) {
      case SUB_JOB:DoSubmitJob();
        break;
      default:break;
    }

    Event *e2p = m_myTimeLine->PopNext();
    if (e2p != currentEvent) {
      cout << "[TGTraceFB::SchedulerAlarmPortal] error: "
           << " e2p != currentEvent " << endl;
    }
    delete currentEvent;
  }

  UpdateAlarm();
}

/* called by simulator */
void
TGTraceFB::NotifyTrafficFinish(double alarmTime,
                               vector<Coflow *> *cfpVp,
                               vector<Flow *> *fpVp) {

  for (vector<Flow *>::iterator fIt = fpVp->begin();
       fIt != fpVp->end(); fIt++) {
    // m_totalFCT += (alarmTime - (*fIt)->GetStartTime());
    if ((*fIt)->GetEndTime() != INVALID_TIME
        && (*fIt)->GetEndTime() >= (*fIt)->GetStartTime()) {
      // a valid start/end time => count FCT
      m_totalFCT += ((*fIt)->GetEndTime() - (*fIt)->GetStartTime());
    } else {
      // sth is wrong.
      cout << " error while calculating flow end time : invalid start/end time"
           << endl;
      exit(-1);
    }
  }

  for (vector<Coflow *>::iterator cfIt = cfpVp->begin();
       cfIt != cfpVp->end(); cfIt++) {
    // count coflow num no matter what.
    m_totalCoflowNum++;

    if ((*cfIt)->IsRejected()) {
      // only count completed coflow.
      // do not consider rejected coflow.
      continue;
    }
    m_total_accepted_coflow_num++;

    double cct = INVALID_TIME;

    if ((*cfIt)->GetEndTime() != INVALID_TIME
        && (*cfIt)->GetEndTime() >= (*cfIt)->GetStartTime()) {
      // a valid end time => count CCT
      cct = (*cfIt)->GetEndTime() - (*cfIt)->GetStartTime();
    } else {
      // sth is wrong.
      cout << " error while calculating coflow end time: invalid start/end time"
           << endl;
      exit(-1);
    }

    if ((*cfIt)->GetDeadlineSec() <= 0) {
      // no deadline.
      m_totalCCT += cct;
    } else {
      // consider deadline.
      if (cct <= ((*cfIt)->GetDeadlineSec() + DEADLINE_ERROR_TOLERANCE)) {
        // met deadline.
        // only count CCT on accepted coflows.
        m_total_met_deadline_num++;
        m_totalCCT += cct;
      }
    }
  }

  // remove the job and coflows
  for (vector<Coflow *>::iterator cfIt = cfpVp->begin();
       cfIt != cfpVp->end(); cfIt++) {
    map<Coflow *, JobDesc *>::iterator jobmapIt = m_coflow2job.find(*cfIt);
    if (jobmapIt == m_coflow2job.end()) {
      cout << "error: can't locate the job from the coflow ptr" << endl;
      cout << "return without clearing up the dead job" << endl;
      return;
    }
    JobDesc *jp2rm = jobmapIt->second;

    bool is_rejected = jp2rm->m_coflow->IsRejected();
    double ddl = jp2rm->m_coflow->GetDeadlineSec();
    double cct = (alarmTime - jp2rm->m_offArrivalTime);
    double lowerbound = jp2rm->m_coflow->GetMaxOptimalWorkSpanInSeconds();
    double elec = (jp2rm->m_coflow->GetMaxPortLoadInBits()
        / (double) DEFAULT_LINK_RATE_BPS);
    double cct_over_lowerbound = cct / lowerbound;
    double cct_over_elec = cct / elec;
    double diff_cct_lowerbound = cct - lowerbound;
    if (diff_cct_lowerbound < 1 / (double) OPTC_BPS) {
      diff_cct_lowerbound = 0.0;
    }

    string title_string, info_string;
    AnalyzeOneCoflow(jp2rm->m_coflow, jp2rm->m_id,
                     jp2rm->m_numSrc, jp2rm->m_numDst,
                     title_string, info_string);
    if (!m_audit_file_title_line) {
      m_audit_file_title_line = true;
      m_jobAuditFile
          << title_string << '\t'
          << "tArr" << '\t'
          << "tFin" << '\t'
          << "rej" << '\t'
          << "ddl" << '\t'
          << "cct" << '\t'            // cct
          << "r-optc" << '\t'         // lowerbound of work span
          << "r-elec" << '\t'         // max aggregated load on a port
          << "diff-optc" << '\t'      // cct - lowerbound of work span
          << endl;
    }
    m_jobAuditFile
        << info_string << '\t'
        << jp2rm->m_offArrivalTime << '\t'
        << alarmTime << '\t'
        << is_rejected << '\t'
        << ddl << '\t'
        << cct << '\t'                      // cct
        << cct_over_lowerbound << '\t'      // lowerbound of work span
        << cct_over_elec << '\t'            // max aggregated load on a port
        << diff_cct_lowerbound << '\t'      // cct - lowerbound of work span
        << endl;

    delete *cfIt;
    delete jp2rm;

    m_coflow2job.erase(jobmapIt);

    vector<JobDesc *>::iterator runJIt = find(m_runningJob.begin(),
                                              m_runningJob.end(),
                                              jp2rm);

    if (runJIt == m_runningJob.end()) {
      cout << "error: the job to submit is not in the "
           << " ready job set!" << endl;
    } else {
      m_runningJob.erase(runJIt);
      m_finishJob.push_back(jp2rm);
    }
  }
}

void
TGTraceFB::NotifySimEnd() {

}

void
TGTraceFB::DoSubmitJob() {

  KickStartReadyJobsAndNotifyScheduler();

  // add next job
  vector<JobDesc *> nextjobs2add = ReadJobs();
  ScheduleToAddJobs(nextjobs2add);
}

void
TGTraceFB::KickStartReadyJobsAndNotifyScheduler() {

  EventSubmitJobDesc *ep = (EventSubmitJobDesc *) m_myTimeLine->PeekNext();
  if (ep->GetEventType() != SUB_JOB) {
    cout << "[TGTraceFB::DoSubmitJob] error: "
         << " the event type is not SUB_JOB!" << endl;
    return;
  }

  vector<JobDesc *>::iterator jIt = find(m_readyJob.begin(),
                                         m_readyJob.end(),
                                         ep->m_jobPtr);

  if (jIt == m_readyJob.end()) {
    cout << "error: the job to submit is not in the "
         << " ready job set!" << endl;
  } else {
    m_readyJob.erase(jIt);
    m_runningJob.push_back(ep->m_jobPtr);
  }

  // dump traffic into network
  vector<Coflow *> *msgCfVp = new vector<Coflow *>();
  msgCfVp->push_back(ep->m_jobPtr->m_coflow);
  MsgEventAddCoflows *msgEp = new MsgEventAddCoflows(m_currentTime, msgCfVp);
  m_simPtr->AddEvent(msgEp);
}

void
TGTraceFB::ScheduleToAddJobs(vector<JobDesc *> &jobs_to_add) {
  if (jobs_to_add.empty()) return;
  for (vector<JobDesc *>::iterator jobIt = jobs_to_add.begin();
       jobIt != jobs_to_add.end(); jobIt++) {
    EventSubmitJobDesc *ep =
        new EventSubmitJobDesc((*jobIt)->m_offArrivalTime, *jobIt);
    m_myTimeLine->AddEvent(ep);
    m_readyJob.push_back(*jobIt);
  }
}

// read next job(s)
// the last job trace line
// should end with return
vector<JobDesc *>
TGTraceFB::ReadJobs() {

  vector<JobDesc *> result;

  string jobLine = "";

  long firstJobTime = -1;

  while (!m_jobTraceFile.eof()) {
    getline(m_jobTraceFile, jobLine);

    if (jobLine.size() <= 0) {
      //cout << "no more jobs are available!" << endl;
      return result;
    }

    vector<string> subFields;
    long numFields = split(jobLine, subFields, '\t');

    if (numFields != 5) {
      cout << "[TGTraceFB::ReadJobs] number of fields illegal!"
           << "Return with job list." << endl;
      return result;
    }

    long jobOffArrivalTimel = stol(subFields[1]);

    if (firstJobTime < 0) {
      firstJobTime = jobOffArrivalTimel;
    }

    if (jobOffArrivalTimel > firstJobTime) {
      // this job should not be read
      //seek back file seeker and return
      m_jobTraceFile.seekg(-(jobLine.length() + 1), m_jobTraceFile.cur);
      return result;
    }

    int jobid = stoi(subFields[0]);
    double jobOffArrivalTime = stod(subFields[1]) / 1000.0;
    int map = stoi(subFields[2]);
    int red = stoi(subFields[3]);

    // perturb if needed.
    // when perturb = false,
    //  if EQUAL_FLOW_TO_SAME_REDUCER = true, all flows to the same reducer
    //     will be the of the same size.
    Coflow *cfp = CreateCoflowPtrFromString(jobOffArrivalTime, jobid,
                                            map, red, subFields[4],
                                            ENABLE_PERTURB_IN_PLAY,
                                            EQUAL_FLOW_TO_SAME_REDUCER);
    cfp->SetJobId(jobid);

    int num_flow = 0;
    if (cfp) {
      num_flow = (int) cfp->GetFlows()->size();
    } else {
      cout << "Error: cfp NULL upon create!" << endl;
    }
    JobDesc *newJobPtr = new JobDesc(jobid,
                                     jobOffArrivalTime,
                                     map, red, num_flow, cfp);
    // add entry into the map
    m_coflow2job.insert(pair<Coflow *, JobDesc *>(cfp, newJobPtr));
    result.push_back(newJobPtr);
  }
  return result;
}

Coflow *
TGTraceFB::CreateCoflowPtrFromString(double time, int coflow_id,
                                     int numMap, int numRed,
                                     string cfInfo,
                                     bool do_perturb, bool avg_size) {

  vector<int> mapLoc = vector<int>(numMap, -1);
  vector<int> redLoc = vector<int>(numRed, -1);
  vector<long> redInput = vector<long>(numRed, (long) 0);

  vector<string> subFields;
  long numFields = split(cfInfo, subFields, '#');

  if (numFields != 2) {
    cout << __func__ << ": number of fields illegal!"
         << "Return with NULL coflow ptr." << endl;
    return NULL;
  }

  vector<string> mapperInfo;
  long numMapFields = split(subFields[0], mapperInfo, ',');
  if (numMapFields != numMap) {
    cout << __func__ << ": numMapFields != numMap "
         << "Return with NULL coflow ptr." << endl;
    return NULL;
  }

  for (int mapCount = 0; mapCount < numMap; mapCount++) {
    mapLoc[mapCount] = stoi(mapperInfo[mapCount]);
  }

  vector<string> reducerInfo;
  long numRedFields = split(subFields[1], reducerInfo, ',');
  if (numRedFields != numRed) {
    cout << __func__ << ": numRedFields != numRed "
         << "Return with NULL coflow ptr." << endl;
    return NULL;
  }

  for (int redCount = 0; redCount < numRed; redCount++) {
    vector<string> reducerDetail;
    long numRedDetailFields = split(reducerInfo[redCount], reducerDetail, ':');
    if (numRedDetailFields != 2) {
      cout << __func__ << ": numRedDetailFields != 2 "
           << "Return withan NULL coflow ptr." << endl;
      return NULL;
    }
    redLoc[redCount] = stoi(reducerDetail[0]);
    redInput[redCount] = 1000000 * stol(reducerDetail[1]);
  }

  vector<Flow *> fpV = vector<Flow *>();
  if (!do_perturb) {
    if (avg_size) {
      CreateFlowsWithEqualSizeToSameReducer(time, numMap, numRed,
                                            mapLoc, redLoc, redInput, fpV);
    } else {
      CreateFlowsWithExactSize(time, numMap, numRed,
                               mapLoc, redLoc, redInput, fpV);
    }
  } else {
    // let us allow some random perturbation.
    unsigned int rand_seed = GetSeedForCoflow(coflow_id);
    CreateFlowsWithSizePerturb(time, numMap, numRed,
                               mapLoc, redLoc, redInput, fpV,
                               5/* hard code of +/-5% */, rand_seed);
  }

  if (fpV.size() <= 0) {
    cout << "error: fpV size <= 0" << endl;
    return NULL;
  }

  int numf = (int) fpV.size();
  Coflow *cfp = new Coflow(time, numf);
  for (vector<Flow *>::iterator fpIt = fpV.begin();
       fpIt != fpV.end(); fpIt++) {
    cfp->AddFlow(*fpIt);
  }
  // initialize the static alpha upon creation.
  cfp->SetStaticAlpha(cfp->CalcAlpha());

  // generate a deadline here if needed!
  if (DEADLINE_MODE) {
    double lb_optc = cfp->GetMaxOptimalWorkSpanInSeconds();
    double lb_elec =
        ((double) cfp->GetLoadOnMaxOptimalWorkSpanInBits()) / ELEC_BPS;
    unsigned int rand_seed = GetSeedForCoflow(coflow_id);
    srand(rand_seed);
    // currently we assume the inflation x = 1;
    double deadline = lb_optc + lb_optc * (((double) (rand() % 100)) / 100.0);
    cfp->SetDeadlineSec(deadline);
    // debug
    cout << " lb_elec " << lb_elec << " lb_optc " << lb_optc << " deadline "
         << deadline << endl;
  }
  return cfp;
}

// use seed_for_seed to initialize the vector of m_seed_for_coflow,
// in length of PERTURB_SEED_NUM.
void TGTraceFB::InitSeedForCoflows(int seed_for_seed) {
  srand(seed_for_seed);
  for (int i = 0; i < PERTURB_SEED_NUM; i++) {
    m_seed_for_coflow.push_back(rand());
    // cout << " seed [" << i << "] = "
    //      << m_seed_for_coflow.back() << endl;
  }
}

int
TGTraceFB::GetSeedForCoflow(int coflow_id) {
  // use jobid to identify a unique seed;
  int seed_idx = (coflow_id + 5122) % PERTURB_SEED_NUM;
  int seed = m_seed_for_coflow[seed_idx];
  // cout << "job id " << jobid << " seed " << seed << endl;
  return seed;
}

// sum( flows to a reducer ) == reducer's input specified in redInput.
void
TGTraceFB::CreateFlowsWithExactSize(double time, int numMap, int numRed,
                                    const vector<int> &mapLoc,
                                    const vector<int> &redLoc,
                                    const vector<long> &redInput,
                                    vector<Flow *> &fpV) {
  for (int redCount = 0; redCount < numRed; redCount++) {

    long redInputTmp = redInput[redCount];

    long avgFlowSize = ceil((double) redInputTmp / (double) numMap);

    for (int mapCount = 0; mapCount < numMap; mapCount++) {

      long flowSize = min(avgFlowSize, redInputTmp);

      int src = mapLoc[mapCount];
      int dst = redLoc[redCount];
      Flow *fp = new Flow(time, src, dst, flowSize);
      fpV.push_back(fp);
      redInputTmp -= flowSize;
    }
  }
}

// divide reducer's input size, specified in redInput, to each of the mapper.
void
TGTraceFB::CreateFlowsWithEqualSizeToSameReducer(double time,
                                                 int numMap,
                                                 int numRed,
                                                 const vector<int> &mapLoc,
                                                 const vector<int> &redLoc,
                                                 const vector<long> &redInput,
                                                 vector<Flow *> &fpV) {
  for (int redCount = 0; redCount < numRed; redCount++) {

    long avgFlowSize = ceil((double) redInput[redCount] / (double) numMap);

    for (int mapCount = 0; mapCount < numMap; mapCount++) {
      long flowSize = avgFlowSize;
      int src = mapLoc[mapCount];
      int dst = redLoc[redCount];
      Flow *fp = new Flow(time, src, dst, flowSize);
      fpV.push_back(fp);
    }
  }
}

void
TGTraceFB::CreateFlowsWithSizePerturb(double time, int numMap, int numRed,
                                      const vector<int> &mapLoc,
                                      const vector<int> &redLoc,
                                      const vector<long> &redInput,
                                      vector<Flow *> &fpV,
                                      int perturb_perc,
                                      unsigned int rand_seed) {
  if (perturb_perc > 100) {
    cout << " Warming : try to perturb the flow sizes with more than 1MB "
         << endl;
  }
  // seed the random generator before we preturb,
  // so that given same traffic trace,
  // we will have the same traffic for different schedulers.
  srand(rand_seed);
  // now we generate traffic.
  for (int redCount = 0; redCount < numRed; redCount++) {
    long redInputTmp = redInput[redCount];

    long avgFlowSize = ceil((double) redInputTmp / (double) numMap);

    for (int mapCount = 0; mapCount < numMap; mapCount++) {

      int perturb_direction = (rand() % 2 == 1) ? 1 : -1;

      // perturb_perc = 5 : (-5%, +5%) flow size , exclusive bound
      double rand_0_to_1 = ((double) rand() / (RAND_MAX));
      double perturb_perc_rand =
          perturb_direction * rand_0_to_1 * (double) perturb_perc / 100.0;
      long flowSize = avgFlowSize * (1 + perturb_perc_rand);

      // only allow flows >= 1MB.
      if (flowSize < 1000000) flowSize = 1000000;

      int src = mapLoc[mapCount];
      int dst = redLoc[redCount];
      Flow *fp = new Flow(time, src, dst, flowSize);
      fpV.push_back(fp);
      redInputTmp -= flowSize;

      // debug
      if (DEBUG_LEVEL >= 10) {
        cout << src << "->" << dst << " flow size after perturb " << flowSize
             << " avgFlowSize " << avgFlowSize
             << " perturb_perc_rand " << perturb_perc_rand << endl;
      }
    }
  }
}

// Analyze the characteristics of a coflow stored in cfp, return the the profile
// info with title_string and info_string.
void
TGTraceFB::AnalyzeOneCoflow(Coflow *cfp, int jobid, int map, int red,
                            string &title_string, string &info_string) {

  // read a coflow. start building up knowledge.
  double total_flow_size_MB = 0.0;
  double min_flow_size_MB = -1.0;
  double max_flow_size_MB = -1.0;

  int flow_num = (int) cfp->GetFlows()->size();
  vector<Flow *> *flowVecPtr = cfp->GetFlows();
  for (vector<Flow *>::iterator fpIt = flowVecPtr->begin();
       fpIt != flowVecPtr->end(); fpIt++) {

    double this_flow_size_MB =
        ((*fpIt)->GetSizeInBit() / (double) 8.0 / 1000000.0);
    total_flow_size_MB += this_flow_size_MB; // each flow >= 1MB
    if (min_flow_size_MB > this_flow_size_MB
        || min_flow_size_MB < 0) {
      min_flow_size_MB = this_flow_size_MB;
    }
    if (max_flow_size_MB < this_flow_size_MB
        || max_flow_size_MB < 0) {
      max_flow_size_MB = this_flow_size_MB;
    }
  }

  double avg_flow_size_MB = total_flow_size_MB / (double) flow_num;
  double lb_optc = cfp->GetMaxOptimalWorkSpanInSeconds();
  long lb_elec = cfp->GetMaxPortLoadInBits();

  string cast_pattern = "";
  if (map == 1 && red == 1) {
    cast_pattern = "1";               // single-flow
  } else if (map > 1 && red == 1) {
    cast_pattern = "m21";            // incast
  } else if (map == 1 && red > 1) {
    cast_pattern = "12m";           // one sender, multiple receiver.
  } else if (map > 1 && red > 1) {
    cast_pattern = "m2m";           // many-to-many
  }

  string bin;
  // double lb_elec_MB = lb_elec / (double) 8.0/ 1000000;
  if (avg_flow_size_MB < 5.0) {
    // short if max flow size < 5MB
    bin += 'S';
  } else {
    bin += 'L';
  }
  if (flow_num <= 50) {
    // narrow if max flow# <= 50
    bin += 'N';
  } else {
    bin += 'W';
  }

  std::stringstream title_stream;
  title_stream
      << "jobid" << '\t'
      << "*cast" << '\t'
      << "bin" << '\t'
      << "map" << '\t'
      << "red" << '\t'
      << "flow#" << '\t'
      << "lb-optc" << '\t'   // work span lower-bound
      << "lb-elec" << '\t'   // max aggregated load on a port
      << "ttl_MB" << '\t'    // total
      << "avg_MB" << '\t' //
      << "min_MB" << '\t'
      << "max_MB";

  std::stringstream info_stream;
  info_stream
      << jobid << '\t'
      << cast_pattern << '\t'
      << bin << '\t'
      << map << '\t'
      << red << '\t'
      << flow_num << '\t'
      << lb_optc << '\t'   // work span lower-bound
      << lb_elec << '\t'   // max aggregated load on a port
      << total_flow_size_MB << '\t'
      << avg_flow_size_MB << '\t'
      << min_flow_size_MB << '\t'
      << max_flow_size_MB;

  title_string = title_stream.str();
  info_string = info_stream.str();
}

///////////////////////////////////////////////////////////////////
///////////// Code for back-to-back coflow              ///////////
///////////// Replay FB trace but ignore real arrival time  ///////
///////////// so that one coflow is served at a time.     /////////
///////////////////////////////////////////////////////////////////

TGFBOnebyOne::TGFBOnebyOne() : TGTraceFB() {
  m_last_coflow_finish_time = 0.0;

  m_comptimeAuditFile.open(COMPTIME_AUDIT_FILE_NAME);
  if (!m_comptimeAuditFile.is_open()) {
    cout << "Error: unable to open file "
         << COMPTIME_AUDIT_FILE_NAME << endl;
    cout << "Now terminate the program" << endl;
    exit(-1);
  }

  m_comptime_audit_title_line = false;
}

TGFBOnebyOne::~TGFBOnebyOne() {
  if (m_comptimeAuditFile.is_open()) {
    m_comptimeAuditFile << "Done" << endl;
    m_comptimeAuditFile.close();
  }
}

void
TGFBOnebyOne::DoSubmitJob() {
  TGTraceFB::KickStartReadyJobsAndNotifyScheduler();
}

// do not apply any destruction here!
void
TGFBOnebyOne::LogComptimeStat(vector<Coflow *> *cfpVp) {
  if (cfpVp->size() != 1) {
    // should only have 1 coflow finished at a time.
    return;
  }
  // obtain comp time stats.
  int jobid = cfpVp->front()->GetJobId();
  CompTimeBreakdown comptime = cfpVp->front()->GetCompTime();

  // look up the job for more information.
  map<Coflow *, JobDesc *>::iterator
      jobmapIt = m_coflow2job.find(cfpVp->front());
  if (jobmapIt == m_coflow2job.end()) {
    cout << "error: can't locate the job from the coflow ptr" << endl;
    cout << "return without clearing up the dead job" << endl;
    return;
  }
  JobDesc *job = jobmapIt->second;
  int num_map = job->m_numSrc;
  int num_red = job->m_numDst;

  if (!m_comptime_audit_title_line) {
    m_comptime_audit_title_line = true;
    m_comptimeAuditFile
        << "jobid" << '\t'
        << "map" << '\t' << "red" << '\t'
        << "sun-ttl" << '\t'
        << "shuffleR" << '\t' << "shuffleS" << '\t'
        << "reserve" << '\t'
        << "solstice-ttl" << '\t'
        << "stuff" << '\t' << "slice" << '\t'
        << "vec-avg" << '\t' << "count" << '\t' << "min" << '\t' << "max"
        << '\t'
        << endl;
  }
  m_comptimeAuditFile
      << jobid << '\t'
      << num_map << '\t'
      << num_red << '\t'
      << comptime.GetSunflowTotalTime() << '\t'
      << comptime.m_time_shuffle_random << '\t'
      << comptime.m_time_sort << '\t'
      << comptime.m_time_reserve << '\t'
      << comptime.GetSolsticeTotalTime() << '\t'
      << comptime.m_time_stuff << '\t'
      << comptime.m_time_slice << '\t'
      << comptime.GetVectorAvgTime() << '\t'
      << comptime.m_time_vector.size() << '\t'
      << comptime.GetVectorMinTime() << '\t'
      << comptime.GetVectorMaxTime() << '\t'
      << endl;
}

void
TGFBOnebyOne::NotifyTrafficFinish(double alarmTime,
                                  vector<Coflow *> *cfpVp,
                                  vector<Flow *> *fpVp) {


  // only output comp time when we play with 1by1.
  // log comptime before calling TGTraceFB::NotifyTrafficFinish()
  // as coflow pointers are destroyed in TGTraceFB::NotifyTrafficFinish().
  LogComptimeStat(cfpVp);

  TGTraceFB::NotifyTrafficFinish(alarmTime, cfpVp, fpVp);

  // add another coflow if any.
  if (!cfpVp->empty()) {
    // we have coflow finish. update last finish time.
    m_last_coflow_finish_time = alarmTime;

    // see whether we have more coflows.
    vector<JobDesc *> nextjob2add = ReadJobs();
    TGTraceFB::ScheduleToAddJobs(nextjob2add);

    // update alarm for scheduler.
    UpdateAlarm();
  }

}

vector<JobDesc *>
TGFBOnebyOne::ReadJobs() {

  vector<JobDesc *> result;

  string jobLine = "";
  getline(m_jobTraceFile, jobLine);

  if (jobLine.size() <= 0) {
    //cout << "no more jobs are available!" << endl;
    return result;
  }

  vector<string> subFields;
  long numFields = split(jobLine, subFields, '\t');

  if (numFields != 5) {
    cout << "[TGTraceFB::ReadJobs] number of fields illegal!"
         << "Return with job list." << endl;
    return result;
  }

  int jobid = stoi(subFields[0]);
  double jobOffArrivalTime = m_last_coflow_finish_time;
  int map = stoi(subFields[2]);
  int red = stoi(subFields[3]);

  // perform perturb if needed.
  // when perturb = false,
  //  if EQUAL_FLOW_TO_SAME_REDUCER = true, all flows to the same reducer
  //     will be the of the same size.
  Coflow *cfp = CreateCoflowPtrFromString(jobOffArrivalTime, jobid,
                                          map, red, subFields[4],
                                          ENABLE_PERTURB_IN_PLAY,
                                          EQUAL_FLOW_TO_SAME_REDUCER);
  cfp->SetJobId(jobid);
  int num_flow = 0;
  if (cfp) {
    num_flow = (int) cfp->GetFlows()->size();
  } else {
    cout << "Error: cfp NULL upon create!" << endl;
  }
  JobDesc *newJobPtr = new JobDesc(jobid,
                                   jobOffArrivalTime,
                                   map, red, num_flow, cfp);
  // add entry into the map
  m_coflow2job.insert(pair<Coflow *, JobDesc *>(cfp, newJobPtr));
  result.push_back(newJobPtr);

  return result;

}

////////////////////////////////////////////////////////////////////////////////
///////////// Code for quick analyze of coflows ONLY                    ////////
////////////////////////////////////////////////////////////////////////////////

void
TGFBAnalyzeOnly::AnalyzeCoflow() {
  string jobLine;
  vector<Coflow *> coflows;

  while (!m_jobTraceFile.eof()) {
    getline(m_jobTraceFile, jobLine);
    if (jobLine.size() <= 0) {
      //cout << "no more jobs are available!" << endl;
      break; // while(!m_jobTraceFile.eof())
    }

    vector<string> subFields;
    long numFields = split(jobLine, subFields, '\t');

    if (numFields != 5) {
      break; // while(!m_jobTraceFile.eof())
    }

    int jobid = stoi(subFields[0]);
    double jobOffArrivalTime = stod(subFields[1]) / 1000.0;
    int map = stoi(subFields[2]);
    int red = stoi(subFields[3]);

    // do perturb if needed.
    // when perturb = false,
    //  if EQUAL_FLOW_TO_SAME_REDUCER = true, all flows to the same reducer
    //     will be the of the same size.
    Coflow *cfp = CreateCoflowPtrFromString(jobOffArrivalTime, jobid,
                                            map, red, subFields[4],
                                            ENABLE_PERTURB_IN_PLAY,
                                            EQUAL_FLOW_TO_SAME_REDUCER);

    if (!cfp) {
      cout << "Error: cfp NULL upon create!" << endl;
      break; // while(!m_jobTraceFile.eof())
    }
    cfp->SetJobId(jobid);
    coflows.push_back(cfp);

    string title_string, info_string;
    TGTraceFB::AnalyzeOneCoflow(cfp,
                                jobid,
                                map,
                                red,
                                title_string,
                                info_string);

    if (!m_analyze_title_line) {
      m_analyze_title_line = true;
      m_jobAuditFile << title_string << endl;
    }

    m_jobAuditFile << info_string << endl;
  }

  // clean up coflow we created for analysis purposes.
  for (vector<Coflow *>::iterator
           coflow = coflows.begin();
       coflow != coflows.end();
       coflow++) {
    delete *coflow;
  }

}


///////////////////////////////////////////////////////////////////////////////////
///////////// Code for quick analyze of network utilization ONLY        ///////////
///////////////////////////////////////////////////////////////////////////////////

void
TGFBUtilizationOnly::LoadAllCoflows(vector<Coflow *> &load_to_me) {
  string jobLine;

  // load all coflows from trace.
  while (!m_jobTraceFile.eof()) {
    getline(m_jobTraceFile, jobLine);
    if (jobLine.size() <= 0) {
      //cout << "no more jobs are available!" << endl;
      break; // while(!m_jobTraceFile.eof())
    }

    vector<string> subFields;
    long numFields = split(jobLine, subFields, '\t');

    if (numFields != 5) {
      break; // while(!m_jobTraceFile.eof())
    }

    int jobid = stoi(subFields[0]);
    double jobOffArrivalTime = stod(subFields[1]) / 1000.0;
    int map = stoi(subFields[2]);
    int red = stoi(subFields[3]);

    // do perturb if needed.
    // when perturb = false,
    //  if EQUAL_FLOW_TO_SAME_REDUCER = true, all flows to the same reducer
    //     will be the of the same size.
    Coflow *cfp = CreateCoflowPtrFromString(jobOffArrivalTime, jobid,
                                            map, red, subFields[4],
                                            ENABLE_PERTURB_IN_PLAY,
                                            EQUAL_FLOW_TO_SAME_REDUCER);

    if (!cfp) {
      cout << "Error: cfp NULL upon create!" << endl;
      break; // while(!m_jobTraceFile.eof())
    }
    cfp->SetJobId(jobid);
    load_to_me.push_back(cfp);

    string title_string, info_string;
    TGTraceFB::AnalyzeOneCoflow(cfp,
                                jobid,
                                map,
                                red,
                                title_string,
                                info_string);

    if (!m_analyze_title_line) {
      m_analyze_title_line = true;
      m_jobAuditFile << title_string << endl;
    }

    m_jobAuditFile << info_string << endl;
  }
}

void
TGFBUtilizationOnly::CleanUpCoflows(vector<Coflow *> &clean_me) {
  // clean up coflow we created for analysis purposes.
  for (vector<Coflow *>::iterator
           coflow = clean_me.begin();
       coflow != clean_me.end();
       coflow++) {
    delete *coflow;
  }
}

void
TGFBUtilizationOnly::AnalyzeUtilization(vector<Coflow *> &coflows) {

  double ONE_GIGA_DOUBLE = 1000000000.0; // 1G
  // Specify the desired link rate under which we calculate the idleness.
  double LINK_RATE_BPS_DOUBLE = ONE_GIGA_DOUBLE;
  // Specify the the factor to inflate the flow sizes, before LoadAllCoflows().
  // Flow sizes are inflated in constructor of Coflow().
  TRAFFIC_SIZE_INFLATE = 1.0;

  LoadAllCoflows(m_coflows);

  double UTILIZATION_SAMPLE_INTERVAL = 0.1; // in second

  priority_queue<UtilizationEvent, vector<UtilizationEvent>,
                 UtilizationEventCompare> util_event_queue;
  double last_coflow_arrival_time = 0;
  // insert all coflow arrival event.
  for (auto coflow : coflows) {
    // skip coflow that only require intra-rack flow(s) and
    // thus does NOT go through the network.
    if (coflow->GetMaxPortLoadInBits() <= 0) {
      continue;//with next coflow
    }
    util_event_queue.push(UtilizationEvent(coflow->GetStartTime(),
                                           coflow->GetJobId(),
                                           UEVENT_ARRIVE));
    last_coflow_arrival_time = coflow->GetStartTime();
  }
  // insert all coflow forget event.
  for (auto coflow : coflows) {
    // skip coflow that only require intra-rack flow(s) and
    // thus does NOT go through the network.
    if (coflow->GetMaxPortLoadInBits() <= 0) {
      continue;//with next coflow
    }
    double when_to_forget = coflow->GetStartTime()
        + coflow->GetMaxPortLoadInBits() / LINK_RATE_BPS_DOUBLE;
    util_event_queue.push(UtilizationEvent(when_to_forget,
                                           coflow->GetJobId(),
                                           UEVENT_FORGET));
  }

  // insert sampling events.
  for (double sample_time = 0.0;
       sample_time < last_coflow_arrival_time;
       sample_time += UTILIZATION_SAMPLE_INTERVAL) {
    util_event_queue.push(UtilizationEvent(sample_time,
                                           -1,
                                           UEVENT_SAMPLE));
  }

  map<int, Coflow *> jobid_to_coflow;
  // build up a map for reference.
  for (auto coflow : coflows) {
    jobid_to_coflow.insert(std::make_pair(coflow->GetJobId(), coflow));
  }

  double last_event_time = 0;
  set<int> jobids_in_memory;
  double weight_sum = 0;
  double weighted_util_sum = 0.0;
  // calculate utilization on the busiest port.
  while (!util_event_queue.empty()) {
    UtilizationEvent event = util_event_queue.top();
    double current_event_time = event.event_time;
    UtilizationEventType event_type = event.event_type;
    int event_jid = event.job_id;

    // cout << "now at " << current_event_time <<
    // " " << event_jid << " " << event_type << endl;
    if (event_type == UEVENT_ARRIVE) {
      jobids_in_memory.insert(event_jid);
      // cout << current_event_time << ": adding " << event_jid << endl;
    } else if (event_type == UEVENT_FORGET) {
      jobids_in_memory.erase(event_jid);
      // cout << current_event_time << ": forget " << event_jid << endl;
      // util_event_queue.pop();
      // continue; // with next event;
    }

    if (event_type != UEVENT_SAMPLE) {
      // not a sample event.
      // no more calculation.
      util_event_queue.pop();
      continue; // with next event;
    }


    // now let's measure the utilization
    // by calculating the coflow demand in memory!
    double util_factor = 0.0;
    {
      // Strategy 2.
      // 0 - 1
      if (!jobids_in_memory.empty()) {
        util_factor = 1.0;
      } else {
        // jobids_in_memory is empty
        util_factor = 0.0;
      }
    }
    if (DEBUG_LEVEL >= 15) {
      cout << " time (" << last_event_time << ", " << current_event_time << ") "
           << current_event_time - last_event_time << " "
           << util_factor << " [";
      for (int jid : jobids_in_memory) {
        cout << jid << ", ";
      }
      cout << "]" << endl;
    }

    double duration = current_event_time - last_event_time;
    weight_sum += duration;
    weighted_util_sum += duration * util_factor;

    if (!m_utilization_title_line) {
      m_utilization_title_line = true;
      m_jobAuditFile
          << "start_time" << '\t'
          << "end_time" << '\t'
          << "duration" << '\t'
          << "utilization" << '\t'
          << endl;
    }

    m_jobAuditFile
        << last_event_time << '\t'
        << current_event_time << '\t'
        << current_event_time - last_event_time << '\t'
        << util_factor << '\t'
        << endl;

    // update and proceed.
    last_event_time = current_event_time;
    util_event_queue.pop();
  }

  cout << " * * * * * * * * " << endl;
  cout << "weighted util " << weighted_util_sum / weight_sum << endl;
  cout << " * * * * * * * * " << endl;
  m_jobAuditFile << "weighted util " << weighted_util_sum / weight_sum << endl;

  // clean up all coflow pointers.
  CleanUpCoflows(m_coflows);
}
