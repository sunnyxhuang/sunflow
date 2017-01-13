//
//  edmond.cc
//
// The code in edmond.{h,cc} comes from online search.

#include "edmond.h"
#include <stdlib.h>
#include <math.h>

/* Global variables */
int *A, *END, *WEIGHT, *NEXTPAIR;
int *MATE, *LINK, *BASE, *NEXTVTX, *LASTVTX, *Y, *NEXT_D, *NEXTEDGE;

int LAST_D, DELTA;

int LASTEDGE[3];

int DUMMYVERTEX, DUMMYEDGE;
int U, V;

int newbase, oldbase, nextbase, stopscan, pairpoint;
int neighbor, nextpoint, newlast;
int newmate, oldmate, oldfirst, firstmate, secondmate;
int f, nextedge, nexte, nextu;

int v, i, e;


/* Assign a pointer link to a vertex.  Edge e joins a vertex in blossom */
/* u to a linked vertex. */

void POINTER(int u, int v, int e) {
  int i, del;

#ifdef EDMOND__DEBUG
  printf("Pointer u,v,e=%d %d %d-%d\n",u,v,END[OPPEDGE(e)],END[e]);
#endif

  LINK[u] = -DUMMYEDGE;
  NEXTVTX[LASTVTX[u]] = DUMMYVERTEX;
  NEXTVTX[LASTVTX[v]] = DUMMYVERTEX;

  if (LASTVTX[u] != u) {
    i = MATE[NEXTVTX[u]];
    del = -SLACK(i) / 2;
  } else del = LAST_D;

  i = u;
  while (i != DUMMYVERTEX) {
    Y[i] += del;
    NEXT_D[i] += del;
    i = NEXTVTX[i];
  }
  if (LINK[v] < 0) {
    LINK[v] = e;
    NEXTPAIR[DUMMYEDGE] = DUMMYEDGE;
    SCAN(v, DELTA);
    return;
  } else {
    LINK[v] = e;
    return;
  }
}

/* Scan each vertex in the blossom whose base is x */

void SCAN(int x, int del) {
  int u, del_e;

#ifdef EDMOND__DEBUG
  printf("Scan del=%d x=%d\n",del,x);
#endif

  newbase = BASE[x];
  stopscan = NEXTVTX[LASTVTX[x]];
  while (x != stopscan) {
    Y[x] += del;
    NEXT_D[x] = LAST_D;
    pairpoint = DUMMYEDGE;
    e = A[x];
    while (e != 0) {
      neighbor = END[e];
      u = BASE[neighbor];
      if (LINK[u] < 0) {
        if (LINK[BMATE (u)] < 0 || LASTVTX[u] != u) {
          del_e = SLACK (e);
          if (NEXT_D[neighbor] > del_e) {
            NEXT_D[neighbor] = del_e;
            NEXTEDGE[neighbor] = e;
          }
        }
      } else if (u != newbase) {
        INSERT_PAIR();
      }
      e = A[e];
    }
    x = NEXTVTX[x];
  }
  NEXTEDGE[newbase] = NEXTPAIR[DUMMYEDGE];
}


/* Expands a blossom.  Fixes up LINK and MATE. */
// useful in Emdond matching
void UNPAIR(int oldbase, int oldmate) {
  int e, newbase, u;

#ifdef EDMOND__DEBUG
  printf("Unpair oldbase, oldmate=%d %d\n",oldbase, oldmate);
#endif

  UNLINK(oldbase);
  newbase = BMATE (oldmate);
  if (newbase != oldbase) {
    LINK[oldbase] = -DUMMYEDGE;
    REMATCH(newbase, MATE[oldbase]);
    if (f == LASTEDGE[1])
      LINK[secondmate] = -LASTEDGE[2];
    else LINK[secondmate] = -LASTEDGE[1];
  }
  e = LINK[oldmate];
  u = BEND (OPPEDGE(e));
  if (u == newbase) {
    POINTER(newbase, oldmate, e);
    return;
  }
  LINK[BMATE (u)] = -e;
  do {
    e = -LINK[u];
    v = BMATE (u);
    POINTER(u, v, -LINK[v]);
    u = BEND (e);
  } while (u != newbase);
  e = OPPEDGE (e);
  POINTER(newbase, oldmate, e);
}


/* changes the matching along an alternating path */
/* firstmate is the first base vertex on the path */
/* edge e is the new matched edge for firstmate   */

// useful in Emdond matching
void REMATCH(int firstmate, int e) {
#ifdef EDMOND__DEBUG
  printf("Rematch firstmate=%d e=%d-%d\n",firstmate, END[OPPEDGE(e)], END[e]);
#endif

  MATE[firstmate] = e;
  nexte = -LINK[firstmate];
  while (nexte != DUMMYEDGE) {
    e = nexte;
    f = OPPEDGE (e);
    firstmate = BEND (e);
    secondmate = BEND (f);
    nexte = -LINK[firstmate];
    LINK[firstmate] = -MATE[secondmate];
    LINK[secondmate] = -MATE[firstmate];
    MATE[firstmate] = f;
    MATE[secondmate] = e;
  }
}


/* unlinks subblossoms in a blossom.  oldbase is the base of the blossom to */
/* be unlinked. */

void UNLINK(int oldbase) {
  int k, j = 1;

#ifdef EDMOND__DEBUG
  printf("Unlink oldbase=%d\n",oldbase);
#endif

  i = NEXTVTX[oldbase];
  newbase = NEXTVTX[oldbase];
  nextbase = NEXTVTX[LASTVTX[newbase]];
  e = LINK[nextbase];
  UL2:
  do {
    nextedge = OPPEDGE (LINK[newbase]);
    for (k = 1; k <= 2; ++k) {
      LINK[newbase] = -LINK[newbase];
      BASE[i] = newbase;
      i = NEXTVTX[i];
      while (i != nextbase) {
        BASE[i] = newbase;
        i = NEXTVTX[i];
      }
      newbase = nextbase;
      nextbase = NEXTVTX[LASTVTX[newbase]];
    }
  } while (LINK[nextbase] == nextedge);
  if (j == 1) {
    LASTEDGE[1] = nextedge;
    j++;
    nextedge = OPPEDGE (e);
    if (LINK[nextbase] == nextedge) {
      goto UL2;
    }
  }
  LASTEDGE[2] = nextedge;

  if (BASE[LASTVTX[oldbase]] == oldbase)
    NEXTVTX[oldbase] = newbase;
  else {
    NEXTVTX[oldbase] = DUMMYVERTEX;
    LASTVTX[oldbase] = oldbase;
  }
}


/* updates numerical bounds for linking paths. */
/* called with LAST_D set to the bound on DELTA for the next search */


// useful in Emdond matching
void SET_BOUNDS() {
  int del;

  for (v = 1; v <= U; ++v) {
    if (LINK[v] < 0 || BASE[v] != v) {
      NEXT_D[v] = LAST_D;
      continue;
    }
    LINK[v] = -LINK[v];
    i = v;
    while (i != DUMMYVERTEX) {
      Y[i] -= DELTA;
      i = NEXTVTX[i];
    }
    f = MATE[v];
    if (f != DUMMYEDGE) {
      i = BEND(f);
      del = SLACK(f);
      while (i != DUMMYVERTEX) {
        Y[i] -= del;
        i = NEXTVTX[i];
      }
    }
    NEXT_D[v] = LAST_D;
  }
}

/* undoes all blossoms to get the final matching */

void UNPAIR_ALL() {
  int u;

  for (v = 1; v <= U; ++v) {
    if (BASE[v] != v || LASTVTX[v] == v)
      continue;
    nextu = v;
    NEXTVTX[LASTVTX[nextu]] = DUMMYVERTEX;
    while (1) {
      u = nextu;
      nextu = NEXTVTX[nextu];
      UNLINK(u);
      if (LASTVTX[u] != u) {
        f = (LASTEDGE[2] == OPPEDGE(e)) ? LASTEDGE[1] : LASTEDGE[2];
        NEXTVTX[LASTVTX[BEND(f)]] = u;
      }
      newbase = BMATE (BMATE(u));
      if (newbase != DUMMYVERTEX && newbase != u) {
        LINK[u] = -DUMMYEDGE;
        REMATCH(newbase, MATE[u]);
      }
      while (LASTVTX[nextu] == nextu && nextu != DUMMYVERTEX)
        nextu = NEXTVTX[nextu];
      if (LASTVTX[nextu] == nextu && nextu == DUMMYVERTEX)
        break;
    }
  }
}


/* set up data structures for weighted match */

/* to add a new type, add new case in SetUp() and a Set_X() routine */

void SetUp(Graph gptr, int type) {
  int i, allocsize;
  Graph g;
  EuclidGraph xy;
  MatrixGraph matg;

  if (type == 1) {
    g = (Graph) gptr;
    U = Degree(g, 0);
    V = NumEdges(g);
  } else if (type == 2) {
    xy = (EuclidGraph) gptr;
    U = xy[0][0];
    V = U * (U - 1) / 2;
  } else if (type == 3) {
    matg = (MatrixGraph) gptr;
    U = matg[0];
    V = U * (U - 1) / 2;
  }

  allocsize = (U + 2 * V + 2) * sizeof(int);
  A = (int *) malloc(allocsize);
  END = (int *) malloc(allocsize);
  WEIGHT = (int *) malloc(allocsize);
  for (i = 0; i < U + 2 * V + 2; i++)
    A[i] = END[i] = WEIGHT[i] = 0;

  if (type == 1) SetStandard(g);
  else if (type == 2) SetEuclid(xy);
  else if (type == 3) SetMatrix(matg);
}

/* set up from Type 1 graph. */

void SetStandard(Graph graph) {
  int elabel, adj_node, i, j;
  int u, v, currentedge;
  Edge edge;

  currentedge = U + 2;
  for (i = 1; i <= U; ++i) {
    edge = FirstEdge(graph, i);
    for (j = 1; j <= Degree(graph, i); ++j) {
      adj_node = EndPoint(edge);
      if (i < adj_node) {
        elabel = ELabel(edge) * 2;
        WEIGHT[currentedge - 1] = WEIGHT[currentedge] = 2 * elabel;
        END[currentedge - 1] = i;
        END[currentedge] = adj_node;
        if (A[i] == 0)
          A[i] = currentedge;
        else {
          u = i;
          v = A[i];
          while (v != 0) {
            if (END[v] > adj_node)
              break;
            u = v;
            v = A[v];
          }
          A[u] = currentedge;
          A[currentedge] = v;
        }
        u = adj_node;
        v = A[u];
        while (v != 0) {
          u = v;
          v = A[v];
        }
        A[u] = currentedge - 1;
        currentedge += 2;
      }
      edge = NextEdge(edge);
    }
  }
}

/* set up from Euclidean graph */

void SetEuclid(EuclidGraph graph) {
  int i, j, currentedge;

  currentedge = U + 2;

  for (i = U; i >= 1; --i)
    for (j = i - 1; j >= 1; --j) {
      WEIGHT[currentedge - 1] = WEIGHT[currentedge]
          = 2 * eucdist2(graph, i, j);
      END[currentedge - 1] = i;
      END[currentedge] = j;
      A[currentedge] = A[i];
      A[i] = currentedge;
      A[currentedge - 1] = A[j];
      A[j] = currentedge - 1;
      currentedge += 2;
    }
}

void SetMatrix(MatrixGraph graph) {
  int i, j, currentedge;

  currentedge = U + 2;

  for (i = U; i >= 1; --i)
    for (j = i - 1; j >= 1; --j) {
      WEIGHT[currentedge - 1] = WEIGHT[currentedge]
          = 2 * graph[j * U + i];
      END[currentedge - 1] = i;
      END[currentedge] = j;
      A[currentedge] = A[i];
      A[i] = currentedge;
      A[currentedge - 1] = A[j];
      A[j] = currentedge - 1;
      currentedge += 2;
    }
}

/* Process an edge linking two linked vertices */
/* Note: global variable v set to the base of one end of the linking edge */
// useful in Emdond matching
void PAIR(int *outcome) {
  int u, w, temp;

#ifdef EDMOND__DEBUG
  printf("Pair v=%d\n",v);
#endif

  e = NEXTEDGE[v];
  while (SLACK(e) != 2 * DELTA)
    e = NEXTPAIR[e];
  w = BEND (e);
  LINK[BMATE (w)] = -e;
  u = BMATE (v);
  while (LINK[u] != -e) {
    LINK[u] = -e;
    if (MATE[w] != DUMMYEDGE) {
      temp = v;
      v = w;
      w = temp;
    }
    v = BLINK (v);
    u = BMATE (v);
  }
  if (u == DUMMYVERTEX && v != w) {
    *outcome = 1;
    return;
  }
  newlast = v;
  newbase = v;
  oldfirst = NEXTVTX[v];
  LINK_PATH(e);
  LINK_PATH(OPPEDGE (e));
  NEXTVTX[newlast] = oldfirst;
  if (LASTVTX[newbase] == newbase)
    LASTVTX[newbase] = newlast;
  NEXTPAIR[DUMMYEDGE] = DUMMYEDGE;
  MERGE_PAIRS(newbase);
  i = NEXTVTX[newbase];
  do {
    MERGE_PAIRS(i);
    i = NEXTVTX[LASTVTX[i]];
    SCAN(i, 2 * DELTA - SLACK(MATE[i]));
    i = NEXTVTX[LASTVTX[i]];
  } while (i != oldfirst);
  *outcome = 0;
  return;
}


/* merges a subblossom's pair list into a new blossom's pair list */
/* v is the base of the previously unlinked subblossom */
/* Note: global variable newbase set to the base of the new blossom */
/* 	called with NEXTPAIR[DUMMYEDGE] pointing to the first edge */
/*		on newbase's pair list */

void MERGE_PAIRS(int v) {
#ifdef EDMOND__DEBUG
  printf("Merge Pairs v=%d\n",v);
#endif

  NEXT_D[v] = LAST_D;
  pairpoint = DUMMYEDGE;
  f = NEXTEDGE[v];
  while (f != DUMMYEDGE) {
    e = f;
    neighbor = END[e];
    f = NEXTPAIR[f];
    if (BASE[neighbor] != newbase)
      INSERT_PAIR();
  }
}


/* links the unlinked vertices in the path P(END[e],newbase) */
/* Note: global variable newbase is set to the base vertex of the new blossom */
/*		newlast is set to the last vertex in newbase's current blossom*/

void LINK_PATH(int e) {
  int u;

#ifdef EDMOND__DEBUG
  printf("Link Path e=%d-%d\n", END[OPPEDGE(e)], END[e]);
#endif

  v = BEND (e);
  while (v != newbase) {
    u = BMATE (v);
    LINK[u] = OPPEDGE (e);
    NEXTVTX[newlast] = v;
    NEXTVTX[LASTVTX[v]] = u;
    newlast = LASTVTX[u];
    i = v;
    BASE[i] = newbase;
    i = NEXTVTX[i];
    while (i != DUMMYVERTEX) {
      BASE[i] = newbase;
      i = NEXTVTX[i];
    }
    e = LINK[v];
    v = BEND (e);
  }
}


/* Update a blossom's pair list. */
/* Note: called with global variable e set to the edge to be inserted. */
/*			neighbor set to the vertex at the end of e */
/*			pairpoint set to the next pair on the pair list */

void INSERT_PAIR() {
  int del_e;

#ifdef EDMOND__DEBUG
  printf("Insert Pair e=%d-%d\n",END[OPPEDGE(e)],END[e]);
#endif

  del_e = SLACK(e) / 2;
  nextpoint = NEXTPAIR[pairpoint];

  while (END[nextpoint] < neighbor) {
    pairpoint = nextpoint;
    nextpoint = NEXTPAIR[nextpoint];
  }
  if (END[nextpoint] == neighbor) {
    if (del_e >= SLACK (nextpoint) / 2)
      return;
    nextpoint = NEXTPAIR[nextpoint];
  }
  NEXTPAIR[pairpoint] = e;
  pairpoint = e;
  NEXTPAIR[e] = nextpoint;
  if (NEXT_D[newbase] > del_e)
    NEXT_D[newbase] = del_e;
}

void AddEdge(Graph g, int n, int m, int label) {
  Edge edge1, edge2;

  edge1 = (Edge) malloc(2 * sizeof(struct edge_ent));
  edge2 = edge1 + 1;

  edge1->label = label;
  edge1->endpoint = m;
  edge1->otheredge = edge2;
  edge1->prevedge = NULL;
  edge1->nextedge = g[n].adj_list;
  if (edge1->nextedge != NULL)
    edge1->nextedge->prevedge = edge1;
  g[n].adj_list = edge1;
  g[n].degree++;

  edge2->label = label;
  edge2->endpoint = n;
  edge2->otheredge = edge1;
  edge2->prevedge = NULL;
  edge2->nextedge = g[m].adj_list;
  if (edge2->nextedge != NULL)
    edge2->nextedge->prevedge = edge2;
  g[m].adj_list = edge2;
  g[m].degree++;
}

Edge FindEdge(Graph graph, int i, int j) {
  Edge edge;

  edge = graph[i].adj_list;
  while (edge != NULL && edge->endpoint != j)
    edge = edge->nextedge;
  if (edge == NULL) return (NULL);
  else return (edge);
}

int RemoveEdge(Graph graph, Edge edge) {
  Edge other;
  int i, j;

  if (edge == NULL) return (0);
  other = edge->otheredge;
  i = other->endpoint;
  j = edge->endpoint;
  graph[i].degree--;
  graph[j].degree--;
  if (edge->prevedge == NULL) {
    graph[i].adj_list = edge->nextedge;
    if (edge->nextedge != NULL)
      edge->nextedge->prevedge = NULL;
  } else if (edge->nextedge == NULL)
    (edge->prevedge)->nextedge = NULL;
  else {
    (edge->nextedge)->prevedge = edge->prevedge;
    (edge->prevedge)->nextedge = edge->nextedge;
  }
  if (other->prevedge == NULL) {
    graph[j].adj_list = other->nextedge;
    if (other->nextedge != NULL)
      other->nextedge->prevedge = NULL;
  } else if (other->nextedge == NULL)
    (other->prevedge)->nextedge = NULL;
  else {
    (other->nextedge)->prevedge = other->prevedge;
    (other->prevedge)->nextedge = other->nextedge;
  }
  free((edge < other) ? edge : other);
  return (1);
}

int NumEdges(Graph g) {
  int i, size, edges;

  edges = 0;
  size = Degree(g, 0);
  for (i = 1; i <= size; i++)
    edges += Degree(g, i);
  edges /= 2;
  return (edges);
}

Graph NewGraph(int size) {
  Graph tmp;
  int i;

  tmp = (Graph) malloc((size + 1) * sizeof(struct node_entry));
  for (i = 1; i <= size; i++) {
    Degree(tmp, i) = 0;
    FirstEdge(tmp, i) = NULL;
    NLabel(tmp, i) = i;
  }
  Degree(tmp, 0) = size;
  return (tmp);
}

EuclidGraph NewEuclid(int size) {
  EuclidGraph xy;

  xy = (EuclidGraph) malloc((size + 1) * 2 * sizeof(int));
  xy[0][0] = size;
  return (xy);
}

MatrixGraph NewMatrix(int size) {
  MatrixGraph graph;
  int i;

  graph = (MatrixGraph) malloc((size * (size + 1) + 1) * sizeof(int));
  graph[0] = size;

  for (i = 1; i <= size; i++)        /* zero the diagonal */
    graph[i * (size + 1)] = 0;

  return (graph);
}

Graph CopyGraph(Graph g) {
  int i, j, size;
  Edge edge;
  Graph cp;

  size = Degree(g, 0);
  cp = NewGraph(size);
  for (i = 1; i <= size; i++) {
    Xcoord(cp, i) = Xcoord(g, i);
    Ycoord(cp, i) = Ycoord(g, i);
    edge = FirstEdge(g, i);
    for (j = 1; j <= Degree(g, i); j++) {
      if (i < EndPoint(edge))
        AddEdge(cp, i, EndPoint(edge), ELabel(edge));
      edge = NextEdge(edge);
    }
  }
  return (cp);
}

//
// useful in Emdond matching
// main entry
int *Weighted_Match(Graph gptr, int type, int maximize) {
  int g, j, w, outcome;
//    int loop=1;

  /* set up internal data structure */
  SetUp(gptr, type);
  Initialize(maximize);

  for (;;) {
    /* printf("Augment #%d\n",loop++); */
    DELTA = 0;
    for (v = 1; v <= U; ++v)
      if (MATE[v] == DUMMYEDGE)
        POINTER(DUMMYVERTEX, v, DUMMYEDGE);
    for (;;) {
      i = 1;
      for (j = 2; j <= U; ++j)
        if (NEXT_D[i] > NEXT_D[j])
          i = j;
      DELTA = NEXT_D[i];
      if (DELTA == LAST_D)
        goto done;
      v = BASE[i];
      if (LINK[v] >= 0) {
        PAIR(&outcome);
        if (outcome == 1)
          break;
      } else {
        w = BMATE (v);
        if (LINK[w] < 0) {
          POINTER(v, w, OPPEDGE(NEXTEDGE[i]));
        } else UNPAIR(v, w);
      }
    }

    LAST_D -= DELTA;
    SET_BOUNDS();
    g = OPPEDGE(e);
    REMATCH(BEND(e), g);
    REMATCH(BEND(g), e);
  }

  done:
  SET_BOUNDS();
  UNPAIR_ALL();
  for (i = 1; i <= U; ++i) {
    MATE[i] = END[MATE[i]];
    if (MATE[i] == DUMMYVERTEX)
      MATE[i] = 0;
  }

  FreeUp();
  return (MATE);
}

// useful in Emdond matching
void Initialize(int maximize) {
  int i, allocsize, max_wt = -MAXWT, min_wt = MAXWT;

  DUMMYVERTEX = U + 1;
  DUMMYEDGE = U + 2 * V + 1;
  END[DUMMYEDGE] = DUMMYVERTEX;

  for (i = U + 2; i <= U + 2 * V; i += 2) {
    if (WEIGHT[i] > max_wt)
      max_wt = WEIGHT[i];
    if (WEIGHT[i] < min_wt)
      min_wt = WEIGHT[i];
  }
  if (!maximize) {
    if (U % 2 != 0) {
      printf("Must have an even number of vertices to do a\n");
      printf("minimum complete matching.\n");
      exit(0);
    }
    max_wt += 2; /* Don't want all zero weight */
    for (i = U + 1; i <= U + 2 * V; i++)
      WEIGHT[i] = max_wt - WEIGHT[i];
    max_wt = max_wt - min_wt;
  }
  LAST_D = max_wt / 2;

  allocsize = (U + 2) * sizeof(int);
  MATE = (int *) malloc(allocsize);
  LINK = (int *) malloc(allocsize);
  BASE = (int *) malloc(allocsize);
  NEXTVTX = (int *) malloc(allocsize);
  LASTVTX = (int *) malloc(allocsize);
  Y = (int *) malloc(allocsize);
  NEXT_D = (int *) malloc(allocsize);
  NEXTEDGE = (int *) malloc(allocsize);
  allocsize = (U + 2 * V + 2) * sizeof(int);
  NEXTPAIR = (int *) malloc(allocsize);

  for (i = 1; i <= U + 1; ++i) {
    MATE[i] = DUMMYEDGE;
    NEXTEDGE[i] = DUMMYEDGE;
    NEXTVTX[i] = 0;
    LINK[i] = -DUMMYEDGE;
    BASE[i] = i;
    LASTVTX[i] = i;
    Y[i] = LAST_D;
    NEXT_D[i] = LAST_D;
  }
}

// useful in Emdond matching
void FreeUp() {
  free(LINK);
  free(BASE);
  free(NEXTVTX);
  free(LASTVTX);
  free(Y);
  free(NEXT_D);
  free(NEXTEDGE);
  free(NEXTPAIR);
  free(A);
  free(END);
  free(WEIGHT);
}

/* Euclidean distance routines */
// kind of useless
int eucdist(EuclidGraph graph, int i, int j) /* Find the distance between two points */
/* 10K x 10K unit square */
{
  int dv, dh;
  register int k, l;

  dv = graph[i][0] - graph[j][0];
  dh = graph[i][1] - graph[j][1];
  k = dv * dv + dh * dh;
  if (k == 0) return (0);
  if (dv < 0) dv = -dv;
  if (dh < 0) dh = -dh;
  l = dv + dh;
  l = (l + k / l) >> 1;
  l = (l + k / l) >> 1;
  l = (l + k / l) >> 1;
  l = (l + k / l) >> 1;
  return ((l * l < k) ? ++l : l);
}

int eucdist2(EuclidGraph graph, int i, int j) /* Find the distance between two points */
/* 1M x 1M unit square */
{
  double dv, dh, d;
  int l;

  dv = (double) graph[i][0] - graph[j][0];
  dh = (double) graph[i][1] - graph[j][1];
  d = sqrt(dv * dv + dh * dh);
  l = (int) d;
  return ((d - l > .000000001) ? l + 1 : l);
}

int eucdistsq(EuclidGraph graph, int i, int j) /* Find the square of the dist between two points */
{
  register int dv, dh;

  dv = graph[i][0] - graph[j][0];
  dh = graph[i][1] - graph[j][1];
  return (dv * dv + dh * dh);
}

///* Graph I/O routines */
//
//Graph ReadGraph (int *size, char file[])
//{	Graph graph;
//    FILE *fp;
//    char c;
//    int edges, degree, vlabel, elabel, adj_node, i, j;
//    int xcoord, ycoord;
//    
//    if (file[0] == '\0') fp = stdin;
//    else fp = fopen(file,"r");
//    if (fp==NULL) {
//        printf("ReadGraph: file %s can't be opened\n",file);
//        exit(0);
//    }
//    fscanf(fp,"%d%d %c",size,&edges,&c);
//    if (c !='U' && c!='u') {
//        printf("ReadGraph: file %s does not contain an undirected graph\n",file);
//        exit(0);
//    }
//    while (getc(fp)!='\n') ;
//    
//    graph = NewGraph(*size);
//    for (i = 1; i <= *size; ++i) {
//        fscanf(fp,"%d%d%d%d",&degree,&vlabel,&xcoord,&ycoord);
//        NLabel(graph,i) = vlabel;
//        Xcoord(graph,i) = xcoord;
//        Ycoord(graph,i) = ycoord;
//        while (getc(fp)!='\n') ;
//        for (j = 1; j <= degree; ++j) {
//            fscanf(fp,"%d%d", &adj_node, &elabel);
//            while (getc(fp)!='\n') ;
//            if (i<adj_node)
//                AddEdge (graph,i,adj_node,elabel);
//        }
//    }
//    fclose(fp);
//    return(graph);
//}
//
//void WriteGraph (Graph graph,char file[])
//{	FILE *fp;
//    int size, i,j,edges;
//    Edge p;
//    
//    if (file[0] == '\0') fp = stdout;
//    else fp = fopen(file,"w");
//    if (fp==NULL) {
//        printf("WriteGraph: file %s can't be opened\n",file);
//        exit(0);
//    }
//    size = Degree(graph,0);
//    edges = NumEdges(graph);
//    fprintf(fp,"%d %d U\n",size,edges);
//    
//    for (i = 1; i <= size; i++) {
//        fprintf(fp,"%d %d %d %d L\n",Degree(graph,i),NLabel(graph,i),
//                Xcoord(graph,i),Ycoord(graph,i));
//        p = FirstEdge(graph,i);
//        for (j = 1; j <= Degree(graph,i); ++j) {
//            fprintf(fp,"%d %d L\n",EndPoint(p),ELabel(p));
//            p = NextEdge(p);
//        }
//    }
//    fclose(fp);
//}
//
//EuclidGraph ReadEuclid(int *size,char file[])
//{	EuclidGraph graph;
//    FILE *fp;
//    char c;
//    int i,xcoord, ycoord;
//    
//    if (file[0]=='\0') fp=stdin;
//    else fp = fopen(file,"r");
//    if (fp==NULL) {
//        printf("ReadEuclid: file %s can't be opened\n",file);
//        exit(0);
//    }
//    fscanf(fp,"%d %c",size,&c);
//    if (c!='E' && c!='e') {
//        printf("ReadEuclid: file %s isn't Euclidean\n",file);
//        exit(0);
//    }
//    while (getc(fp)!='\n');
//    graph = NewEuclid(*size);
//    
//    for (i=1; i<=*size; ++i) {
//        fscanf(fp,"%d%d",&xcoord,&ycoord);
//        while (getc(fp)!='\n') ;
//        graph[i][0] = xcoord;
//        graph[i][1] = ycoord;
//    }
//    fclose(fp);
//    return (graph);
//}
//
//void WriteEuclid(EuclidGraph graph, char file[])
//{	FILE *fp;
//    int size, i;
//    
//    if (file[0] == '\0') fp = stdout;
//    else fp = fopen(file,"w");
//    if (fp==NULL) {
//        printf("WriteEuclid: file %s can't be opened\n",file);
//        exit(0);
//    }
//    size = graph[0][0];
//    fprintf(fp,"%d E\n",size);
//    for (i = 1; i <= size; i++)
//        fprintf(fp,"%d %d\n",graph[i][0],graph[i][1]);
//    
//    fclose(fp);
//}
//
//MatrixGraph ReadMatrix(int *size, char file[])
//{	MatrixGraph graph;
//    FILE *fp;
//    char c;
//    int i,j,k;
//    
//    if (file[0]=='\0') fp=stdin;
//    else fp = fopen(file,"r");
//    if (fp==NULL) {
//        printf("ReadMatrix: file %s can't be opened\n",file);
//        exit(0);
//    }
//    fscanf(fp,"%d %c",size,&c);
//    if (c!='M' && c!='m') {
//        printf("ReadMatrix: file %s isn't a distance matrix\n",file);
//        exit(0);
//    }
//    while (getc(fp)!='\n');
//    graph = NewMatrix(*size);
//    
//    for (i=1; i<*size; i++) {
//        for (j=i+1; j<=*size; j++) {
//            fscanf(fp,"%d",&k);
//            graph[i*(*size)+j] = graph[j*(*size)+i] = k;
//        }
//        while (getc(fp)!='\n');
//    }
//    fclose(fp);
//    return(graph);
//}
//
//void WriteMatrix(MatrixGraph graph, char file[])
//{	FILE *fp;
//    int size, i, j;
//    
//    if (file[0] == '\0') fp = stdout;
//    else fp = fopen(file,"w");
//    if (fp==NULL) {
//        printf("WriteMatrix: file %s can't be opened\n",file);
//        exit(0);
//    }
//    size = graph[0];
//    fprintf(fp,"%d M\n",size);
//    for (i = 1; i < size; i++) {
//        for (j=i+1; j<=size; j++)
//            fprintf(fp,"%d ",graph[i*size+j]);
//        fprintf(fp,"\n");
//    }
//    fclose(fp);
//}



