//
//  main.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 9/20/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
// 

#include <iostream>

#include "events.h"

using namespace std;

int main(int argc, const char *argv[]) {

  // *- "sunflow"  : Sunflow circuit scheduler
  //    note that when reconf = 0, shut down resource leakage check
  //                            or change to a larger tolerance.

  // *- "varysImpl" : varys implemention _similar_ to that on GitHub
  //                  = selfish coflow + flow based Work Conservation
  //                  Tweaked by Sunny because the performance of the original
  //                  implementation is too bad when compared with Sunflow.
  //                  After tweaked, total/avg CCT performance is improved and
  //                  becomes comparable with Sunflow.
  //  - "varysImpl201606" : varys implemention as seen on GitHub in 2016 June.
  //                        Total/avg CCT performance is much worse than Sunflow.
  //  - "aaloImpl" : aalo implemented as seen in GitHub

  // *- "Solstice" : Conext'15. Big-slice-decomp with threshold of halved upper
  //                 bound.

  //  - "BvN" : optimal under circuit reconfiguration delay = 0.
  //  - "TMS" : optc fair share + TMS (improved with Solstice's quick_stuff)
  //  - "Edmond" : optc fair share + Edmond

  string schedulerName = "sunflow";

  TRAFFIC_SIZE_INFLATE = 1;
  DEBUG_LEVEL = 0;

  // "fb" | "fb1by1" | "utilizationONLY" | "analyzeONLY"
  string trafficProducerName = "fb1by1";

  for (int i = 1; i < argc; i = i + 2) {
    /* We will iterate over argv[] to get the parameters stored inside.
     * Note that we're starting on 1 because we don't need to know the
     * path of the program, which is stored in argv[0] */
    if (i + 1 != argc) { // Check that we haven't finished parsing already
      string strFlag = string(argv[i]);
      if (strFlag == "-elec") {
        string content(argv[i + 1]);
        ELEC_BPS = stol(content);
      } else if (strFlag == "-optc") {
        string content(argv[i + 1]);
        OPTC_BPS = stol(content);
      } else if (strFlag == "-d") {
        string content(argv[i + 1]);
        CIRCUIT_RECONF_DELAY = stod(content);
      } else if (strFlag == "-inflate") {
        string content(argv[i + 1]);
        TRAFFIC_SIZE_INFLATE = stod(content);
      } else if (strFlag == "-traffic") {
        trafficProducerName = string(argv[i + 1]);
      } else if (strFlag == "-ftrace") {
        TRAFFIC_TRACE_FILE_NAME = string(argv[i + 1]);
      } else if (strFlag == "-faudit") {
        TRAFFIC_AUDIT_FILE_NAME = string(argv[i + 1]);
      } else if (strFlag == "-caudit") {
        COMPTIME_AUDIT_FILE_NAME = string(argv[i + 1]);
      } else if (strFlag == "-paudit") {
        PORT_AUDIT_FILE_NAME = string(argv[i + 1]);
      } else if (strFlag == "-naudit") {
        NET_TX_AUDIT_FILE_NAME = string(argv[i + 1]);
      } else if (strFlag == "-ciraudit") {
        CIRCUIT_AUDIT_FILE_NAME = string(argv[i + 1]);
      } else if (strFlag == "-s") {
        schedulerName = string(argv[i + 1]);
      } else if (strFlag == "-zc") {
        string content(argv[i + 1]);
        if (content == "true") {
          ZERO_COMP_TIME = true;
        } else {
          ZERO_COMP_TIME = false;
        }
      } else if (strFlag == "-sunflow_shuffle_random") {
        string content(argv[i + 1]);
        if (content == "true") {
          SUNFLOW_SHUFFLE_RANDOM = true;
        } else {
          SUNFLOW_SHUFFLE_RANDOM = false;
        }
      } else if (strFlag == "-sunflow_sort") {
        string content(argv[i + 1]);
        if (content == "true") {
          SUNFLOW_SORT_DEMAND = true;
        } else {
          SUNFLOW_SORT_DEMAND = false;
        }
      } else {
        cout << "invalid arguments " << strFlag << " \n";
        exit(0);
      }
    }
  }

  // DEFAULT_LINK_RATE_BPS for bottle-neck calculation
  DEFAULT_LINK_RATE_BPS = OPTC_BPS > ELEC_BPS ? OPTC_BPS : ELEC_BPS;

  Simulator ximulator;
  ximulator.InstallScheduler(schedulerName);
  ximulator.InstallTrafficGen(trafficProducerName);

  cout << "schedulerName = " << schedulerName << endl;
  cout << "trafficProducer = " << trafficProducerName << endl;
  cout << "TRAFFIC_SIZE_INFLATE = " << TRAFFIC_SIZE_INFLATE << endl;
  cout << "ELEC_BPS = " << ELEC_BPS << endl;
  cout << "OPTC_BPS = " << OPTC_BPS << endl;
  cout << "CIRCUIT_RECONF_DELAY = " << CIRCUIT_RECONF_DELAY << endl;
  cout << "ZERO_COMP_TIME = " << std::boolalpha << ZERO_COMP_TIME << endl;
  cout << "NUM_RACK = " << NUM_RACK << " *  NUM_LINK_PER_RACK = "
       << NUM_LINK_PER_RACK << endl;
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << "ENABLE_PERTURB_IN_PLAY = " << std::boolalpha
       << ENABLE_PERTURB_IN_PLAY << endl;
  cout << "LOG_COMP_STAT = " << std::boolalpha << LOG_COMP_STAT << endl;
  cout << " *  *  " << endl;
  cout << "AUDIT_PORT_ALL_SRC = " << std::boolalpha << PORT_AUDIT_ALL_SRC
       << endl;
  cout << "AUDIT_PORT_LAST_GONE = "
       << std::boolalpha << PORT_AUDIT_LAST_GONE_SRC_DST << endl;
  cout << "AUDIT_PORT_BOTTLENECK = "
       << std::boolalpha << PORT_AUDIT_BOTTLENECK << endl;
  cout << " *  *  " << endl;
  int file_name_cutoff = (int) TRAFFIC_TRACE_FILE_NAME.size() - 30;
  if (file_name_cutoff < 0) file_name_cutoff = 0;
  cout << "TRAFFIC_TRACE_FILE_NAME = "
       << TRAFFIC_TRACE_FILE_NAME.substr(file_name_cutoff) << endl;

  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << " (1MB size/optc_rate/switching) = "
       << 1000000.0 * 8 / OPTC_BPS / CIRCUIT_RECONF_DELAY << endl;
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << "Sunflow only: SUNFLOW_SHUFFLE_RANDOM = " << std::boolalpha
       << SUNFLOW_SHUFFLE_RANDOM << endl;
  cout << "Sunflow only: SUNFLOW_SORT_DEMAND = " << std::boolalpha
       << SUNFLOW_SORT_DEMAND << endl;
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << "TMS only: TMS_CYCLE_LENGTH = " << TMS_CYCLE_LENGTH << endl;
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";
  cout << " *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * \n";

  ximulator.Run();

  return 0;
}


