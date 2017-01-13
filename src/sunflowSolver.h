//
//  sunflowSolver.h
//  Ximulator
//
//  Created by Xin Sunny Huang on 2/11/16.
//  Copyright © 2016 Xin Sunny Huang. All rights reserved.
//

#ifndef sunflowSolver_h
#define sunflowSolver_h

#include <map>
#include <set>
#include <vector>

using namespace std;

// helper class
typedef enum Port_Event_Type {
  kIDLE,
  kBUSY,
  kUNKNOWN
} PortEventType;

typedef enum Port_Type {
  kSRC,
  kDST
} PortType;

struct PortEvent {
 public:
  PortEvent() : m_event_type(kUNKNOWN),
                m_start_time(-1.0), m_end_time(-1.0) {}

  PortEvent(PortEventType event_type, double start_time, double end_time,
            PortType self_type, int self_port,
            PortType peer_port_type, int peer_port) :
      m_event_type(event_type), m_start_time(start_time), m_end_time(end_time),
      m_port_type(self_type), m_port(self_port),
      m_peer_port_type(peer_port_type), m_peer_port(peer_port) {}

  PortEventType m_event_type;
  double m_start_time;
  double m_end_time;
  PortType m_port_type;
  int m_port;
  // required if event_type = kBUSY;
  PortType m_peer_port_type;
  int m_peer_port;
};

class PortEventQueue {
 public:
  PortEventQueue() {}
  ~PortEventQueue() { m_port_event_queue.clear(); }

  void Push(PortEvent port_event);
  PortEvent Pop();
  double Peek();
  void Clear() { m_port_event_queue.clear(); }
  bool Empty() { return m_port_event_queue.empty(); }
  unsigned long Size() { return m_port_event_queue.size(); }
 private:
  vector<PortEvent> m_port_event_queue;
};

class SunflowSolver {
 public:
  SunflowSolver();
  virtual ~SunflowSolver() {};

  void ReserveCircuits(vector<map<pair<int, int>, long>> &layered_demand);

  bool HasNext() { return !m_port_reserve.Empty(); };
  double GetNextCircuitsAndExpireTime(map<pair<int, int>,
                                          double> &circuits_to_expire_time);
  double GetMaxReservationTime() { return m_max_reserv_time; }

  // get comptime stats.
  double GetCompTimeReserve();
 private:

  void ReserveCircuitsInternal(vector<map<pair<int, int>,
                                          long>> &layered_demand);

  void ReserveCircuitsFirstCoflow(vector<pair<pair<int, int>,
                                              long>> &demand_vector,
                                  double &max_reserv_time);

  void ReserveCircuitsOtherCoflow(vector<pair<pair<int, int>,
                                              long>> &demand_vector,
                                  double max_reserv_time);

  //used to perform connectivity test. No-op by default.
  void GetSrcDstSocketCandidates(
      int src, int dst,
      const set<int> &busy_src_sockets,
      const set<int> &busy_dst_sockets,
      vector<pair<int, int>> &src_dst_socket_candidates);

  virtual void ConnectivityRest() {}
  virtual bool ConnectivityTest(int src_socket, int dst_socket,
                                const set<int> &busy_src_sockets,
                                const set<int> &busy_dst_sockets) { return true; }
  virtual void ConnectivityAddEdge(int src, int dst_socket) {}
  virtual void ConnectivityRemoveSocketSrc(int src_socket) {}
  virtual void ConnectivityRemoveSocketDst(int dst_socket) {}

  // used by sunflow to shuffle demand entries.
  void ShuffleDemand(vector<pair<pair<int, int>, long>> &demand_vector);
  bool IsSingleSrcOrSingleDst(const vector<pair<pair<int, int>,
                                                long>> &demand_vector);

  void BuildPortFreeUpsFromReservations(double max_reserv_time);
  void InitializeBusyCircuitsAndFreeUntil(set<int> &busy_src,
                                          set<int> &busy_dst,
                                          map<int, double> &src_free_until,
                                          map<int, double> &dst_free_until);

  // reservations on src w/ dst peer info
  PortEventQueue m_port_reserve;

  // freeups on dst and src
  PortEventQueue m_port_freeups;

  // resulting parameters.
  double m_max_reserv_time;

  // setup parameters.
  // transcribed from global settings.
  long m_optc_bps;
  double m_reconf_delay;
  bool m_enable_shuffle_random;
  bool m_enable_shuffle_sorted_demand;

  // computation stat.
  vector<double> m_shuffle_time;
  vector<double> m_reserve_time;

  // debug
  void PrintReservations(const PortEventQueue &);

 protected:
  // setup parameters.
  // transcribed from global settings.
  int m_debug_level;
  int m_num_rack;
  int m_num_link_per_rack;
};

#endif /* sunflowSolver_h */
