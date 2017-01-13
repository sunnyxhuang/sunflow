//
//  global.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 10/30/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//

#include "global.h"
#include "util.h"

// Some useful constant/scaler
// const double ONE_GIGA_DOUBLE = 1000000000.0;
const long ONE_GIGA_LONG = 1000000000;
// const long TEN_GIGA_LONG = 10000000000;

long DEFAULT_LINK_RATE_BPS = ONE_GIGA_LONG / NUM_LINK_PER_RACK;
long ELEC_BPS = ONE_GIGA_LONG / NUM_LINK_PER_RACK; //1G = 1000000000
long OPTC_BPS = ONE_GIGA_LONG / NUM_LINK_PER_RACK; //10G = 10000000000

double CIRCUIT_RECONF_DELAY = 0.01000000; //0.01000000 ＝ 10ms // 0.00001000 = 10us

// used by the scheduler to hack the simulation if needed.
bool SPEEDUP_1BY1_TRAFFIC_IN_SCHEDULER = false;
bool DEADLINE_MODE = false;
double DEADLINE_ERROR_TOLERANCE = 0.0001;

bool ENABLE_LEAK_CHECK = false;
bool LEAK_CHECK_EXIT = false;

bool ENABLE_PERTURB_IN_PLAY = true;
// valid if ENABLE_PERTURB_IN_PLAY = false;
// all flows to the same reducer will be equally
// distributed for all mappers.
bool EQUAL_FLOW_TO_SAME_REDUCER = false;

// used to initialized end time.
double INVALID_TIME = -1.0;

// valid for sunflow scheduler only. 
// if true, shuffle demand entries (after shuffle by decomp) by random before scheduling.
bool SUNFLOW_SHUFFLE_RANDOM = false;
// if true, sort demand entries (after shuffle) before scheduling.
bool SUNFLOW_SORT_DEMAND = false;

// if true, log down comp time stat.
bool LOG_COMP_STAT = false;

// for tms scheduling only.
// when larger than 0, only the largest N slots will be applied.
int TMS_LARGEST_N_SLOTS_ONLY = 10;
// tms schedule a number of slot into one scheduling cycle.
double TMS_CYCLE_LENGTH = 0.2;
// for sinkhorn's algorithm
long SH_EPSLON = 1;
double SH_ERROR_BOUND = 0.001 / 10000;
long SH_MAX_ITR_NUM = 1000000000;

// edmond schedule one slot per scheduling cycle.
double EDMOND_CYCLE_LENGTH = 0.02;

// default num of racks used in tms to determine
// the bound for random selected rack to fill demand
int NUM_RACK = 150;
//  NUM_LINK_PER_RACK circuits are connected for a rack;
//  rack_number = circuit_number mod NUM_RACK
int NUM_LINK_PER_RACK = 1;

// used by traffic generator.
// which indicates max number of coflows supported per run.
// used to initialized the random seed for each coflow.
int PERTURB_SEED_NUM = 600;

// all flow byte sizes are multiplied by TRAFFIC_SIZE_INFLATE.
double TRAFFIC_SIZE_INFLATE = 1;

// TODO: replace ${YOUR_DIR_TO_SUNFLOW} with your directory to the sunflow
// package, e.g. "/Users/yourid/Documents/research/sunflow"
// string SRC_DIR = "${YOUR_DIR_TO_SUNFLOW}/";
string SRC_DIR = "../"; // this only works if you run your binary under src/
string TRAFFIC_TRACE_FILE_NAME = SRC_DIR + "trace/fbtrace-1hr.txt";
string RESULTS_DIR = "results/";
string TRAFFIC_AUDIT_FILE_NAME = SRC_DIR + RESULTS_DIR + "audit_traffic.txt";
string COMPTIME_AUDIT_FILE_NAME = SRC_DIR + RESULTS_DIR +  "comp_time.txt";
string PORT_AUDIT_FILE_NAME = SRC_DIR + RESULTS_DIR + "audit_port.txt";
string NET_TX_AUDIT_FILE_NAME = SRC_DIR + RESULTS_DIR + "audit_net_tx.txt";
string CIRCUIT_AUDIT_FILE_NAME = SRC_DIR + RESULTS_DIR + "audit_circuits.txt";

// selections of port auditing.
bool PORT_AUDIT_LAST_GONE_SRC_DST = false;
bool PORT_AUDIT_BOTTLENECK = false;
bool PORT_AUDIT_ALL_SRC = false;

bool FILL_DEMAND = true;
bool SORT_SLOTS_REAL_DEMAND = true;
bool ZERO_COMP_TIME = true;
bool VALIDATE_DSM = true;

// for aalo.
int AALO_Q_NUM = 10;
double AALO_INIT_Q_HEIGHT = 10.0 * 1000000; // 10MB
double AALO_Q_HEIGHT_MULTI = 10.0;

// a flag to indicated this coflow has inf large alpha (expected cct).
double CF_DEAD_ALPHA_SIGN = -1.0;
// consider online alpha to be infinite if longer than this cutoff,
// so that the coflow may be skipped and enter the so-called "fair-share" trick.
double ONLINE_ALPHA_CUTOFF = 10000000; // 1000000000

// by default 0 - no debug string
// more debug string with higher DEBUG_LEVEL
int DEBUG_LEVEL = 0;

// if true, ignore switching circuit set warning; false by default.
bool IGNORE_WARNING_SWITCHING_CIRCUIT_SET = false;
