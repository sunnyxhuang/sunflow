//
//  hopcroft.cc
//
// The code in hopcroft.{h,cc} comes from online search.
//

#include "hopcroft.h" 
#include <string.h>
#include <iostream>

int n1, n2, edges, last[MAXN1], prev[MAXM], head[MAXM];
int matching[MAXN2], dist[MAXN1], Q[MAXN1];
bool used[MAXN1], vis[MAXN1];

void hopcroftInit(int _n1, int _n2) {
  n1 = _n1;
  n2 = _n2;
  edges = 0;
  memset(last, -1, sizeof(int) * n1);
}

void hopcroftAddEdge(int u, int v) {
  head[edges] = v;
  prev[edges] = last[u];
  last[u] = edges++;
}

void bfs() {
  memset(dist, -1, sizeof(int) * n1);
  int sizeQ = 0;
  int u, i, e;
  for (u = 0; u < n1; ++u) {
    if (!used[u]) {
      Q[sizeQ++] = u;
      dist[u] = 0;
    }
  }
  for (i = 0; i < sizeQ; i++) {
    int u1 = Q[i];
    for (e = last[u1]; e >= 0; e = prev[e]) {
      int u2 = matching[head[e]];
      if (u2 >= 0 && dist[u2] < 0) {
        dist[u2] = dist[u1] + 1;
        Q[sizeQ++] = u2;
      }
    }
  }
}

bool dfs(int u1) {
  vis[u1] = true;
  int e;
  for (e = last[u1]; e >= 0; e = prev[e]) {
    int v = head[e];
    int u2 = matching[v];
    if (u2 < 0 || (!vis[u2] && dist[u2] == dist[u1] + 1 && dfs(u2))) {
      matching[v] = u1;
      used[u1] = true;
      return true;
    }
  }
  return false;
}

int hopcroft(int *match_result) {
  memset(used, 0, sizeof(int) * n1);
  memset(matching, -1, sizeof(int) * n2);
  int res, u;
  for (res = 0;;) {
    bfs();
    memset(vis, 0, sizeof(int) * n1);
    int f = 0;
    for (u = 0; u < n1; ++u)
      if (!used[u] && dfs(u))
        ++f;
    if (!f) {
      memcpy(match_result, matching, sizeof(int) * n2);
      return res;
    }
    res += f;
  }
}
