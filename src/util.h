//
//  util.h
//  Ximulator
//
//  Created by Xin Sunny Huang on 9/21/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//

#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <set>
#include <vector>
#include <string>
#include <time.h>
using namespace std;

int const FLOAT_TIME_WIDTH = 10;

unsigned long split(const string &txt, vector<string> &strs, char ch);

double secondPass(struct timeval end_time, struct timeval start_time);

bool CompByFlowValue(std::pair<pair<int, int>, long> const &a,
                     std::pair<pair<int, int>, long> const &b);
template<typename K, typename V>
V MapWithDef(std::map<K, V> &m, const K &key, const V &defval) {
  typename std::map<K, V>::const_iterator it = m.find(key);
  if (it == m.end()) {
    //m[key] = defval;
    m.insert(pair<K, V>(key, defval));
    return defval;
  } else {
    return it->second;
  }
}

template<typename K, typename V>
V FindWithDef(const std::map<K, V> &m, const K &key, const V &defval) {
  typename std::map<K, V>::const_iterator it = m.find(key);
  if (it == m.end()) {
    //m[key] = defval;
    //m.insert(pair<K, V>(key, defval));
    return defval;
  } else {
    return it->second;
  }
}

template<typename K, typename V>
V MapWithInc(std::map<K, V> &m, const K &key, const V &inc) {
  typename std::map<K, V>::const_iterator it = m.find(key);
  if (it == m.end()) {
    //m[key] = inc;
    m.insert(pair<K, V>(key, inc));
  } else {
    m[key] += inc;
  }
  return m[key];
}

template<typename K, typename V>
V MinMap(std::map<K, V> &m) {
  if (m.empty()) {
    return NULL;
  }
  typename std::map<K, V>::iterator it;
  V minVal = m.begin()->second;
  for (it = m.begin(); it != m.end(); it++) {
    if (minVal > it->second) {
      minVal = it->second;
    }
  }
  return minVal;
}

template<typename K, typename V>
K Key2MinPositiveMap(std::map<K, V> &m) {
  if (m.empty()) {
    return NULL;
  }
  typename std::map<K, V>::iterator it;
  V minVal = m.begin()->second;
  K minKey = m.begin()->first;
  for (it = m.begin(); it != m.end(); it++) {
    if ((minVal <= 0 || minVal > it->second)
        && it->second > 0) {
      // positive
      minKey = it->first;
      minVal = it->second;
    }
  }
  return minKey;
}

template<typename K, typename V>
K Key2MaxPositiveMap(std::map<K, V> &m) {
  if (m.empty()) {
    return 0;
  }
  typename std::map<K, V>::iterator it;
  V maxVal = m.begin()->second;
  K maxKey = m.begin()->first;
  for (it = m.begin(); it != m.end(); it++) {
    if ((maxVal <= 0 || maxVal < it->second)
        && it->second > 0) {
      // positive
      maxKey = it->first;
      maxVal = it->second;
    }
  }
  return maxKey;
}

//return the max value in a map
template<typename K, typename V>
V MaxMap(std::map<K, V> &m) {
  if (m.empty()) {
    return 0.0;
  }
  typename std::map<K, V>::iterator it;
  V maxVal = m.begin()->second;
  for (it = m.begin(); it != m.end(); it++) {
    if (maxVal < it->second) {
      maxVal = it->second;
    }
  }
  return maxVal;
}

template<typename K>
bool ContainsKey(const std::set<K> &s, const K &key) {
  typename std::set<K>::const_iterator it = s.find(key);
  return it != s.end();
}

template<typename K, typename V>
bool ContainsKey(const std::map<K, V> &m, const K &key) {
  typename std::map<K, V>::const_iterator it = m.find(key);
  return it != m.end();
}

#endif /*UTIL_H*/
