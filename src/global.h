//
//  global.h
//  Ximulator
//
//  Created by Xin Sunny Huang on 10/27/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//

#ifndef Ximulator_global_h
#define Ximulator_global_h

#include <iostream>

using namespace std;

extern long DEFAULT_LINK_RATE_BPS;
extern long ELEC_BPS;
extern long OPTC_BPS;

extern double CIRCUIT_RECONF_DELAY;

extern bool SPEEDUP_1BY1_TRAFFIC_IN_SCHEDULER;
extern bool DEADLINE_MODE;
extern double DEADLINE_ERROR_TOLERANCE;

extern bool ENABLE_LEAK_CHECK;
extern bool LEAK_CHECK_EXIT;
extern bool ENABLE_PERTURB_IN_PLAY;
extern bool EQUAL_FLOW_TO_SAME_REDUCER;

extern double INVALID_TIME;

extern bool SUNFLOW_SHUFFLE_RANDOM;
extern bool SUNFLOW_SORT_DEMAND;

extern bool LOG_COMP_STAT;

extern int TMS_LARGEST_N_SLOTS_ONLY;
extern double TMS_CYCLE_LENGTH;
extern long SH_EPSLON;
extern double SH_ERROR_BOUND;
extern long SH_MAX_ITR_NUM;

extern double EDMOND_CYCLE_LENGTH;

extern int NUM_RACK;
extern int NUM_LINK_PER_RACK;

extern int PERTURB_SEED_NUM;

extern string TRAFFIC_TRACE_FILE_NAME;
extern string TRAFFIC_AUDIT_FILE_NAME;
extern string COMPTIME_AUDIT_FILE_NAME;

// port "utilization"
extern string PORT_AUDIT_FILE_NAME;
extern bool PORT_AUDIT_LAST_GONE_SRC_DST;
extern bool PORT_AUDIT_BOTTLENECK;
extern bool PORT_AUDIT_ALL_SRC;
// for network usage measurement
extern string NET_TX_AUDIT_FILE_NAME;

extern bool FILL_DEMAND;
extern bool SORT_SLOTS_REAL_DEMAND;
extern bool ZERO_COMP_TIME;
extern bool VALIDATE_DSM;

extern double TRAFFIC_SIZE_INFLATE;

// for aalo.
extern int AALO_Q_NUM;
extern double AALO_INIT_Q_HEIGHT;
extern double AALO_Q_HEIGHT_MULTI;

// a flag to indicated this coflow has inf large alpha.
extern double CF_DEAD_ALPHA_SIGN;
extern double ONLINE_ALPHA_CUTOFF;

// for circuit auditing.
extern string CIRCUIT_AUDIT_FILE_NAME;

extern int DEBUG_LEVEL;
extern bool IGNORE_WARNING_SWITCHING_CIRCUIT_SET;

typedef enum Event_Type {
  NO_EVENT,
  SUB_JOB,                /*for TrafficGen*/
  WORKER_FINISH,
  FLOW_FINISH,            /*for Scheduler */
  COFLOW_ARRIVE,
  FLOW_ARRIVE,
  RESCHEDULE,
  ORDER_CIRCUIT,
  APPLY_CIRCUIT,
  APPLY_NEW_SCHEDULE,
  SCHEDULE_END,
  MSG_ADD_FLOWS,          /*for Simulator */
  MSG_ADD_COFLOWS,
  ALARM_TRAFFIC,
  MSG_TRAFFIC_FINISH,
  ALARM_SCHEDULER,
} EventType;


#endif
