//
//  sunflowSolver.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 2/11/16.
//  Copyright © 2016 Xin Sunny Huang. All rights reserved.
//

#include "global.h"
#include "util.h"
#include "sunflowSolver.h"

#include <algorithm>
#include <iostream>
#include <sys/time.h>

SunflowSolver::SunflowSolver() {
  m_max_reserv_time = -1;
  m_debug_level = DEBUG_LEVEL;
  m_optc_bps = OPTC_BPS;
  m_reconf_delay = CIRCUIT_RECONF_DELAY;
  m_enable_shuffle_random = SUNFLOW_SHUFFLE_RANDOM;
  m_enable_shuffle_sorted_demand = SUNFLOW_SORT_DEMAND;
  m_num_rack = NUM_RACK;
  m_num_link_per_rack = NUM_LINK_PER_RACK;
}

double
SunflowSolver::GetNextCircuitsAndExpireTime(map<pair<int, int>, double> &
circuits_to_expire_time) {
  circuits_to_expire_time.clear();

  if (m_port_reserve.Empty()) {
    return 0.0;
  }

  double first_reservation_time = -1;

  while (!m_port_reserve.Empty()) {

    double this_reserv_time = m_port_reserve.Peek();

    if (this_reserv_time <= first_reservation_time
        || first_reservation_time < 0) {
      first_reservation_time = this_reserv_time;
    } else {
      // first_reservation_time < this_reserv_time
      break; // no more qualifying reservations.
    }

    PortEvent this_reservation = m_port_reserve.Pop();
    int src = this_reservation.m_port;
    int dst = this_reservation.m_peer_port;
    double expire_time = this_reservation.m_end_time;
    circuits_to_expire_time.insert(std::make_pair(std::make_pair(src, dst),
                                                  expire_time));
  }
  return first_reservation_time;
}

void
SunflowSolver::ReserveCircuits(vector<map<pair<int, int>,
                                          long>> &layered_demand) {

  ReserveCircuitsInternal(layered_demand);

  //debug
  PrintReservations(m_port_reserve);
}

void
SunflowSolver::ReserveCircuitsInternal(vector<map<pair<int, int>,
                                                  long>> &layered_demand) {
  m_port_reserve.Clear();
  m_port_freeups.Clear();

  if (layered_demand.empty()) return;

  double local_max_reserv_time = -1.0;
  for (vector<map<pair<int, int>, long>>::iterator
           one_coflow = layered_demand.begin();
       one_coflow != layered_demand.end();
       one_coflow++) {
    // for each coflow

    vector<pair<pair<int, int>, long>> demand_vector(one_coflow->begin(),
                                                     one_coflow->end());

    /* ----------------  shuffle ------------------------------- */
    struct timeval shuffle_start_time;
    gettimeofday(&shuffle_start_time, NULL);

    ShuffleDemand(demand_vector);

    struct timeval shuffle_end_time;
    gettimeofday(&shuffle_end_time, NULL);
    m_shuffle_time.push_back(secondPass(shuffle_end_time, shuffle_start_time));

    /* ----------------  reserve ------------------------------- */
    struct timeval reserve_start_time;
    gettimeofday(&reserve_start_time, NULL);

    if (local_max_reserv_time < 0) {
      // take care of the first coflow
      ReserveCircuitsFirstCoflow(demand_vector,
                                 local_max_reserv_time);
    } else {
      // take care of the rest coflows
      ReserveCircuitsOtherCoflow(demand_vector,
                                 local_max_reserv_time);
    }

    struct timeval reserve_end_time;
    gettimeofday(&reserve_end_time, NULL);
    m_reserve_time.push_back(secondPass(reserve_end_time, reserve_start_time));
  }
  m_max_reserv_time = local_max_reserv_time;
}

void
SunflowSolver::GetSrcDstSocketCandidates(
    int src,
    int dst,
    const set<int> &busy_src_sockets,
    const set<int> &busy_dst_sockets,
    vector<pair<int, int> > &src_dst_socket_candidates) {

  src_dst_socket_candidates.clear();

  vector<int> src_socket_candidates; // circuit number at the src port
  vector<int> dst_socket_candidates; // circuit number at the dst port
  for (int i = 0; i < m_num_link_per_rack; i++) {
    if (!ContainsKey(busy_src_sockets, src + i * m_num_rack)) {
      src_socket_candidates.push_back(src + i * m_num_rack);
    }
  }
  for (int i = 0; i < m_num_link_per_rack; i++) {
    if (!ContainsKey(busy_dst_sockets, dst + i * m_num_rack)) {
      dst_socket_candidates.push_back(dst + i * m_num_rack);
    }
  }

  if (src_socket_candidates.empty() || dst_socket_candidates.empty()) {
    return;
  }

  for (vector<int>::const_iterator
           this_src_socket = src_socket_candidates.begin();
       this_src_socket != src_socket_candidates.end();
       this_src_socket++) {

    for (vector<int>::const_iterator
             this_dst_socket = dst_socket_candidates.begin();
         this_dst_socket != dst_socket_candidates.end();
         this_dst_socket++) {

      if (ConnectivityTest(*this_src_socket, *this_dst_socket,
                           busy_src_sockets, busy_dst_sockets)) {
        // setting up this circuit agrees with the connectivity requirements.
        src_dst_socket_candidates.push_back(std::make_pair(*this_src_socket,
                                                          *this_dst_socket));
      }
    }
  }
}

void
SunflowSolver::ReserveCircuitsFirstCoflow(vector<pair<pair<int, int>,
                                                      long>> &demand_vector,
                                          double &max_reserv_time) {

  double expect_finish_time = 0.0;

  set<int> busy_src_sockets;
  set<int> busy_dst_sockets;
  ConnectivityRest();

  double circuit_available_time = 0.0;

  while (!demand_vector.empty()) {

    for (vector<pair<pair<int, int>, long>>::iterator
             demand_iter = demand_vector.begin();
         demand_iter != demand_vector.end();) {
      int src = demand_iter->first.first;
      int dst = demand_iter->first.second;

      // select candidates.
      vector<pair<int, int>> src_dst_socket_candidates;
      GetSrcDstSocketCandidates(src, dst,
                                busy_src_sockets, busy_dst_sockets,
                                src_dst_socket_candidates);

      if (!src_dst_socket_candidates.empty()) {
        int src_socket = src_dst_socket_candidates.front().first;
        int dst_socket = src_dst_socket_candidates.front().second;

        // compatible circuit
        busy_src_sockets.insert(src_socket);
        busy_dst_sockets.insert(dst_socket);
        ConnectivityAddEdge(src_socket, dst_socket);

        double circuit_duration = m_reconf_delay
            + demand_iter->second / (double) m_optc_bps;
        double circuit_start_time = circuit_available_time;
        double circuit_freeup_time = circuit_available_time + circuit_duration;

        // add reservations
        PortEvent reserve_src_dst = PortEvent(kBUSY,
                                              circuit_start_time,
                                              circuit_freeup_time,
                                              kSRC, src_socket,
                                              kDST, dst_socket);

        m_port_reserve.Push(reserve_src_dst);

        if (m_debug_level >= 15
            || (m_debug_level >= 1 && (src == 1 || dst == 64))) {
          cout << "FIRST coflow reserves a circuit "
               << src_socket << " -> " << dst_socket
               << " (" << src << " -> " << dst << ") "
               << " during (" << circuit_start_time << ", "
               << circuit_freeup_time << ")"
               << endl;
        }

        // add freeups
        PortEvent freeup_src_dst = PortEvent(kIDLE,
                                             circuit_freeup_time,
                                             -1,
                                             kSRC,
                                             src_socket,
                                             kDST,
                                             dst_socket);
        m_port_freeups.Push(freeup_src_dst);

        // remove demand
        demand_iter = demand_vector.erase(demand_iter);

        // update finish time
        if (expect_finish_time < circuit_freeup_time) {
          expect_finish_time = circuit_freeup_time;
        }

      } else {
        demand_iter++;
      }
    }

    if (!demand_vector.empty()) {
      // pop an circuit to free
      PortEvent freeup_src_dst = m_port_freeups.Pop();
      circuit_available_time = freeup_src_dst.m_start_time;
      int src_socket = freeup_src_dst.m_port;
      int dst_socket = freeup_src_dst.m_peer_port;
      busy_src_sockets.erase(src_socket);
      busy_dst_sockets.erase(dst_socket);
      ConnectivityRemoveSocketSrc(src_socket);
      ConnectivityRemoveSocketDst(dst_socket);
    }

  } // while(!demand_vector.empty())
  max_reserv_time = expect_finish_time;

  // one_coflow_compstat.m_time_reserve = time_reserve;
}

//return true if there is only one src or dst.
bool SunflowSolver::
IsSingleSrcOrSingleDst(const vector<pair<pair<int, int>,
                                         long>> &demand_vector) {
  set<int> all_src;
  set<int> all_dst;
  for (vector<pair<pair<int, int>, long>>::const_iterator
           demand = demand_vector.begin();
       demand != demand_vector.end();
       demand++) {
    all_src.insert(demand->first.first);
    all_dst.insert(demand->first.second);
  }
  if (all_src.size() <= 1 || all_dst.size() <= 1) {
    return true;
  }
  return false;
}

void
SunflowSolver::ShuffleDemand(
    vector<pair<pair<int, int>, long>> &demand_vector) {

  bool no_need_for_shuffle = IsSingleSrcOrSingleDst(demand_vector);
  if (no_need_for_shuffle) {
    // there is no need for shuffle or sort
    // if (single dst) or (single src)
    return;
  }

  if (m_enable_shuffle_random) {
    // do valid work.
    std::random_shuffle(demand_vector.begin(), demand_vector.end());
    return;
  }

  if (m_enable_shuffle_sorted_demand) {
    std::stable_sort(demand_vector.begin(),
                     demand_vector.end(),
                     CompByFlowValue);
    return;
  }
}

void
SunflowSolver::ReserveCircuitsOtherCoflow(
    vector<pair<pair<int, int>, long>> &demand_vector,
    double max_reserv_time) {

  BuildPortFreeUpsFromReservations(max_reserv_time);

  set<int> busy_src_sockets;
  set<int> busy_dst_sockets;
  map<int, double> src_free_until;
  map<int, double> dst_free_until;

  double circuit_available_time = 0.0;

  ConnectivityRest();
  InitializeBusyCircuitsAndFreeUntil(busy_src_sockets, busy_dst_sockets,
                                     src_free_until, dst_free_until);

  while (!demand_vector.empty()) {

    for (vector<pair<pair<int, int>, long>>::iterator
             demand_iter = demand_vector.begin();
         demand_iter != demand_vector.end();) {

      int src = demand_iter->first.first;
      int dst = demand_iter->first.second;

      // debug
      // cout << "see demand " << src << "->" << dst << endl;

      // debug
      if (m_debug_level >= 1 && src == 1 && dst == 64) {
        cout << "see demand " << src << "->" << dst << endl;
      }

      vector<pair<int, int>> src_dst_socket_candidates;
      GetSrcDstSocketCandidates(src, dst,
                                busy_src_sockets, busy_dst_sockets,
                                src_dst_socket_candidates);

      if (src_dst_socket_candidates.empty()) {
        // no circuit can meet the connectivity requirement.
        demand_iter++;
        continue;
      }

      int src_circuit = src_dst_socket_candidates.front().first;
      double src_feasible_until = MapWithDef(src_free_until,
                                             src_circuit,
                                             max_reserv_time);
      int dst_circuit = src_dst_socket_candidates.front().second;
      double dst_feasible_until = MapWithDef(dst_free_until,
                                             dst_circuit,
                                             max_reserv_time);
      double feasible_until = min(src_feasible_until, dst_feasible_until);
      for (vector<pair<int, int>>::iterator
               src_dst_pair = src_dst_socket_candidates.begin() + 1;
           src_dst_pair != src_dst_socket_candidates.end();
           src_dst_pair++) {
        double this_src_free_until = MapWithDef(src_free_until,
                                                src_dst_pair->first,
                                                max_reserv_time);
        double this_dst_free_until = MapWithDef(dst_free_until,
                                                src_dst_pair->second,
                                                max_reserv_time);
        double this_circuit_free_until = min(this_src_free_until,
                                             this_dst_free_until);
        if (this_circuit_free_until > feasible_until) {
          feasible_until = this_circuit_free_until;
          src_circuit = src_dst_pair->first;
          dst_circuit = src_dst_pair->second;
        }
      }

      double min_required_until =
          circuit_available_time + m_reconf_delay + 1 / (double) m_optc_bps;

      if (feasible_until < min_required_until) {
        // skip small slot.
        demand_iter++;
        continue; // for each demand
      }

      double desired_until
          = circuit_available_time + m_reconf_delay
              + demand_iter->second / (double) m_optc_bps;

      double valid_until = min(feasible_until, desired_until);

      if (valid_until < min_required_until) {
        // demand is close to 0.
        // cout << "see error of circuit_duration" << endl;
        demand_iter++;
        continue; // for each demand
      }

      // see opportunity : compatible circuit & slot
      busy_src_sockets.insert(src_circuit);
      busy_dst_sockets.insert(dst_circuit);
      ConnectivityAddEdge(src_circuit, dst_circuit);

      double circuit_start_time = circuit_available_time;
      double circuit_freeup_time = valid_until;

      // add reservations
      PortEvent reserve_src_dst =
          PortEvent(kBUSY, circuit_start_time, circuit_freeup_time,
                    kSRC, src_circuit, kDST, dst_circuit);
      m_port_reserve.Push(reserve_src_dst);

      if (m_debug_level >= 15 ||
          (m_debug_level >= 1 && (src == 1 || dst == 64))) {
        cout << "OTHER Coflow reserves a circuit "
             << src_circuit << " -> " << dst_circuit
             << " (" << src << " -> " << dst << ") "
             << " during (" << circuit_start_time << ", " << circuit_freeup_time
             << ")"
             << endl;
      }


      // add freeups
      if (valid_until < src_feasible_until) {
        // see free slot opportunity for src.
        PortEvent freeup_src =
            PortEvent(kIDLE, circuit_freeup_time, src_feasible_until,
                      kSRC, src_circuit, kSRC, src_circuit);
        m_port_freeups.Push(freeup_src);
      }
      if (valid_until < dst_feasible_until) {
        // see free slot opportunity for dst.
        PortEvent freeup_dst =
            PortEvent(kIDLE, circuit_freeup_time, dst_feasible_until,
                      kDST, dst_circuit, kDST, dst_circuit);

        m_port_freeups.Push(freeup_dst);
      }

      // remove demand
      if ((desired_until - valid_until) * m_optc_bps <= 0) {
        demand_iter = demand_vector.erase(demand_iter);
      } else {
        demand_iter->second = (desired_until - valid_until) * m_optc_bps;
        demand_iter++;
      }

    } // for each demand

    if (m_port_freeups.Empty()) {
      // no more free lunch.
      break; // while(demand_vector.empty())
    }

    if (!demand_vector.empty()) {
      bool is_last_event_reservation = true;
      while (!m_port_freeups.Empty() &&
          (is_last_event_reservation
              || m_port_freeups.Peek() == circuit_available_time)) {
        // pop until a port is free
        // and pop all events that happens at the same time as
        // the first freeup
        PortEvent port_event = m_port_freeups.Pop();

        int self_port = port_event.m_port;
        int peer_port = port_event.m_peer_port;

        if (m_debug_level >= 1) {
          if (self_port == 1 || peer_port == 0)
            cout << "popping  "
                 << port_event.m_port << "->" << port_event.m_peer_port
                 << " event type " << port_event.m_event_type
                 << " (" << port_event.m_start_time << ", "
                 << port_event.m_end_time << ")"
                 << " freeups left: " << m_port_freeups.Size()
                 << endl;
        }

        if (port_event.m_event_type == kIDLE) {
          // this is a freeup event.
          circuit_available_time = port_event.m_start_time;
          double free_until = port_event.m_end_time;
          // freeup event is only valid for a port
          if (port_event.m_port_type == kSRC) {
            busy_src_sockets.erase(self_port);
            src_free_until[self_port] = free_until;
            ConnectivityRemoveSocketSrc(self_port);
          } else if (port_event.m_port_type == kDST) {
            busy_dst_sockets.erase(self_port);
            dst_free_until[self_port] = free_until;
            ConnectivityRemoveSocketDst(self_port);
          }
          is_last_event_reservation = false;
        } else if (port_event.m_event_type == kBUSY) {
          // this is a reservation event.
          // reservation event is valid for two ports.
          busy_src_sockets.insert(self_port);
          busy_dst_sockets.insert(peer_port);
          ConnectivityAddEdge(self_port, peer_port);
          // will continue to pop.
        } else {
          cout << "internal consistency error while popping port events"
               << endl;
        }
      }
    }
  }
}

void
SunflowSolver::InitializeBusyCircuitsAndFreeUntil(
    set<int> &busy_src,
    set<int> &busy_dst,
    map<int, double> &src_free_until,
    map<int, double> &dst_free_until) {
  while (!m_port_freeups.Empty()
      && m_port_freeups.Peek() == 0) {
    // no freeups at time 0
    PortEvent reservation = m_port_freeups.Pop();

    //debug
    if (m_debug_level >= 15 && reservation.m_port == 70) {
      cout << __func__ << " see port event " << reservation.m_port
           << " -> " << reservation.m_peer_port
           << " during (" << reservation.m_start_time << ", "
           << reservation.m_end_time << ")"
           << endl;
    }

    if (reservation.m_event_type == kBUSY) {

      if (reservation.m_port_type != kSRC
          || reservation.m_peer_port_type != kDST) {
        cout << "internal consistency error while initializing busy src/dst"
             << endl;
        continue;
      }

      // initialize busy circuit.
      busy_src.insert(reservation.m_port);
      busy_dst.insert(reservation.m_peer_port);
      ConnectivityAddEdge(reservation.m_port, reservation.m_peer_port);

    } else if (reservation.m_event_type == kIDLE) {
      // initialize free until.
      double free_until = reservation.m_end_time;
      if (reservation.m_port_type == kSRC) {
        src_free_until[reservation.m_port] = free_until;
      } else if (reservation.m_port_type == kDST) {
        dst_free_until[reservation.m_port] = free_until;
      } else {
        cout << "internal consistency error while initializing free until"
             << endl;
        continue;
      }
    } else {
      cout << __func__ << ": internal consistency error" << endl;
    }

  }
}

// build free up events from existing reservations.
void
SunflowSolver::BuildPortFreeUpsFromReservations(double max_reserv_time) {
  if (DEBUG_LEVEL >= 5) {
    cout << "building freeups from " << m_port_reserve.Size() << " reservations"
         << endl;
  }
  m_port_freeups.Clear();

  if (m_port_reserve.Empty()) return;

  // map src/dst -> (free_start_time, free_end_time)
  map<int, vector<pair<double, double>>> src_free_slots;
  map<int, vector<pair<double, double>>> dst_free_slots;

  PortEventQueue existing_reservations = m_port_reserve;

  while (!existing_reservations.Empty()) {

    PortEvent one_reservation = existing_reservations.Pop();

    // add reservations events into the freeup queue.
    m_port_freeups.Push(one_reservation);

    int src = one_reservation.m_port;
    int dst = one_reservation.m_peer_port;
    double reservation_start_time = one_reservation.m_start_time;
    double reservation_end_time = one_reservation.m_end_time;

    if (m_debug_level >= 15) {
      cout << "see reservation "
           << src << "->" << dst
           << " (" << reservation_start_time << ", " << reservation_end_time
           << ")"
           << endl;

    }

    if (src_free_slots[src].empty()) {
      if (reservation_start_time > 0) {
        src_free_slots[src].push_back(std::make_pair(0,
                                                     reservation_start_time));
      }
      // free slot's end time is initialized as the max_reserv_time.
      src_free_slots[src].push_back(std::make_pair(reservation_end_time,
                                                   max_reserv_time));
    } else if (src_free_slots[src].back().first == reservation_start_time) {
      // vector not empty
      // last free slot start at the beginning of a new reservation
      // no free opportunities
      src_free_slots[src].back().first = reservation_end_time;
    } else if (src_free_slots[src].back().first < reservation_start_time) {
      // we see a real free slot
      // free slot's end time is updated here.
      src_free_slots[src].back().second = reservation_start_time;
      src_free_slots[src].push_back(std::make_pair(reservation_end_time,
                                                   max_reserv_time));
    } else {
      cout << "src " << src << " frees at " << src_free_slots[src].back().first
           << "s but a reservation starts at " << reservation_start_time << "s"
           << endl;
      cout << "internal consistency error when constructing src_free_slots "
           << endl;
    }

    if (dst_free_slots[dst].empty()) {
      if (reservation_start_time > 0) {
        dst_free_slots[dst].push_back(std::make_pair(0,
                                                     reservation_start_time));
      }
      // free slot's end time is initialized as the max_reserv_time.
      dst_free_slots[dst].push_back(std::make_pair(reservation_end_time,
                                                   max_reserv_time));
    } else if (dst_free_slots[dst].back().first == reservation_start_time) {
      // vector not empty
      // last free slot start at the beginning of a new reservation
      // no free opportunities
      dst_free_slots[dst].back().first = reservation_end_time;
    } else if (dst_free_slots[dst].back().first < reservation_start_time) {
      // we see a real free slot
      // free slot's end time is updated here.
      dst_free_slots[dst].back().second = reservation_start_time;
      dst_free_slots[dst].push_back(std::make_pair(reservation_end_time,
                                                   max_reserv_time));
    } else {
      cout << "dst " << dst << "frees at " << dst_free_slots[dst].back().first
           << "s but a reservation starts at " << reservation_start_time << "s"
           << endl;
      cout << "internal consistency error when constructing dst_free_slots "
           << endl;
    }

  };


  // now we construct port freeup events according to src_free_slots and dst_free_slots;
  for (map<int, vector<pair<double, double>>>::const_iterator
           src_freeup_iter = src_free_slots.begin();
       src_freeup_iter != src_free_slots.end();
       src_freeup_iter++) {
    int src = src_freeup_iter->first;
    const vector<pair<double, double>>
        &src_free_slot_list = src_freeup_iter->second;
    for (vector<pair<double, double>>::const_iterator
             one_free_slot = src_free_slot_list.begin();
         one_free_slot != src_free_slot_list.end();
         one_free_slot++) {
      if (one_free_slot->second - one_free_slot->first > 0) {
        // see free slot opportunity.
        double freeup_start_time = one_free_slot->first;
        double freeup_end_time = one_free_slot->second;
        //introduce a bug here
        // freeup_end_time = one_free_slot->second - one_free_slot->first;
        m_port_freeups.Push(PortEvent(kIDLE, freeup_start_time, freeup_end_time,
                                      kSRC, src, kSRC, src));
        if (m_debug_level >= 1 && src == 1) {
          // only work for src
          //debug
          cout << "adding a freeup for " << src
               << " during (" << one_free_slot->first << ", "
               << one_free_slot->second << ")"
               << endl;
        }

      }
    }

  }

  for (map<int, vector<pair<double, double>>>::const_iterator
           dst_freeup_iter = dst_free_slots.begin();
       dst_freeup_iter != dst_free_slots.end();
       dst_freeup_iter++) {
    int dst = dst_freeup_iter->first;
    const vector<pair<double, double>>
        &dst_free_slot_list = dst_freeup_iter->second;
    for (vector<pair<double, double>>::const_iterator
             one_free_slot = dst_free_slot_list.begin();
         one_free_slot != dst_free_slot_list.end();
         one_free_slot++) {
      if (one_free_slot->second - one_free_slot->first > 0) {
        // see free slot opportunity.
        double freeup_start_time = one_free_slot->first;
        double freeup_end_time = one_free_slot->second;
        m_port_freeups.Push(PortEvent(kIDLE, freeup_start_time, freeup_end_time,
                                      kDST, dst, kDST, dst));
        if (m_debug_level >= 1 && dst == 64) {
          // debug
          // only work for dst
          cout << "adding a freeup for " << dst
               << " during (" << one_free_slot->first << ", "
               << one_free_slot->second << ")"
               << endl;
        }
      }
    }
  }
}

double
SunflowSolver::GetCompTimeReserve() {
  double sum = 0.0;
  for (vector<double>::const_iterator
           t = m_reserve_time.begin();
       t != m_reserve_time.end();
       t++) {
    sum += *t;
  }
  return sum; // total
}

void
SunflowSolver::PrintReservations(const PortEventQueue &reservations_to_print) {
  if (m_debug_level < 30) {
    return;
  }

  PortEventQueue existing_reservations = reservations_to_print;

  cout << "Print Reservation " << endl;

  while (!existing_reservations.Empty()) {

    PortEvent one_reservation = existing_reservations.Pop();

    int src = one_reservation.m_port;
    int dst = one_reservation.m_peer_port;
    double reservation_start_time = one_reservation.m_start_time;
    double reservation_end_time = one_reservation.m_end_time;

    cout << src << "->" << dst
         << " (" << reservation_start_time
         << ", " << reservation_end_time << ")"
         << endl;
  }
}

// helper class.
void
PortEventQueue::Push(PortEvent port_event) {
  double add_time = port_event.m_start_time;
  for (vector<PortEvent>::iterator it = m_port_event_queue.begin();
       it != m_port_event_queue.end();
       it++) {
    if (add_time < it->m_start_time) {
      m_port_event_queue.insert(it, port_event);
      return;
    }
  }
  m_port_event_queue.push_back(port_event);
}

double
PortEventQueue::Peek() {
  if (m_port_event_queue.empty()) return -1;
  return m_port_event_queue.front().m_start_time;
}

PortEvent
PortEventQueue::Pop() {
  if (m_port_event_queue.empty()) return PortEvent();
  PortEvent to_return = *m_port_event_queue.begin();
  m_port_event_queue.erase(m_port_event_queue.begin());
  return to_return;
}