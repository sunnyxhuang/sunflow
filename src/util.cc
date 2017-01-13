//
//  util.cc
//  Ximulator
//
//  Created by Xin Sunny Huang on 10/19/14.
//  Copyright (c) 2014 Xin Sunny Huang. All rights reserved.
//


#include <string>
#include <ctgmath>

#include "util.h"

using namespace std;

unsigned long split(const string &txt, vector<string> &strs, char ch) {
  unsigned long pos = txt.find(ch);
  unsigned long initialPos = 0;
  strs.clear();

  // Decompose statement
  while (pos != std::string::npos) {
    strs.push_back(txt.substr(initialPos, pos - initialPos));
    initialPos = pos + 1;

    pos = txt.find(ch, initialPos);
  }

  // Add the last one
  unsigned long lastPos = pos < txt.size() ? pos : txt.size();
  strs.push_back(txt.substr(initialPos, lastPos - initialPos));

  return strs.size();
}

// comparator that sorts items large -> small by value(second)
// maintain order if values equal.
bool CompByFlowValue(std::pair<pair<int, int>, long> const &a,
                     std::pair<pair<int, int>, long> const &b) {
  return a.second > b.second;
};

double secondPass(struct timeval end_time, struct timeval start_time) {
  return (double) (end_time.tv_sec - start_time.tv_sec)
      + ((double) (end_time.tv_usec - start_time.tv_usec)) / (double) 1000000;
}
 