//
//  edmond.h
//
// The code in edmond.{h,cc} comes from online search.

#ifndef __Ximulator__edmond__
#define __Ximulator__edmond__

#include <stdio.h>

#define    MAXWT  100000000
#define INF    1000000000 //1000000000 // 2,147,483,647?

/* the number of the blossom entered by edge e */
#define BEND(e) (BASE[END[e]])

/* the blossom matched with v's blossom */
#define BMATE(v) (BASE[END[MATE[v]]])

/* the blossom entered by the edge that links v's blossom */
#define BLINK(v) (BASE[END[LINK[v]]])

/* the edge e with it's direction reversed */
#define OPPEDGE(e) (((e - U) % 2 == 0) ? (e - 1) : (e + 1))

/* the slack of edge e */
#define SLACK(e) (Y[END[e]] + Y[END[OPPEDGE(e)]] - WEIGHT[e])

/* Global variables */
extern int *A, *END, *WEIGHT, *NEXTPAIR;
extern int *MATE, *LINK, *BASE, *NEXTVTX, *LASTVTX, *Y, *NEXT_D, *NEXTEDGE;

extern int LAST_D, DELTA;

extern int LASTEDGE[3];

extern int DUMMYVERTEX, DUMMYEDGE;
extern int U, V;

extern int newbase, oldbase, nextbase, stopscan, pairpoint;
extern int neighbor, nextpoint, newlast;
extern int newmate, oldmate, oldfirst, firstmate, secondmate;
extern int f, nextedge, nexte, nextu;

extern int v, i, e;

struct node_entry {
  int degree;
  int label;
  int x;
  int y;
  struct edge_ent *adj_list;
};
typedef struct node_entry *Graph;

struct edge_ent {
  int endpoint;
  int label;
  int label2;
  struct edge_ent *nextedge;
  struct edge_ent *prevedge;
  struct edge_ent *otheredge;
};
typedef struct edge_ent *Edge;

#define Degree(graph, n)    (graph[n].degree)
#define NLabel(graph, n)    (graph[n].label)
#define Xcoord(graph, n)    (graph[n].x)
#define Ycoord(graph, n)    (graph[n].y)
#define FirstEdge(graph, n) (graph[n].adj_list)

#define EndPoint(e) (e->endpoint)
#define ELabel(e)   (e->label)
#define ELabel2(e)  (e->label2)
#define Other(e)    (e->otheredge)
#define NextEdge(e) (e->nextedge)

/* Euclidean graph type */
typedef int (*EuclidGraph)[2];

/* Distance matrix graph type */
typedef int *MatrixGraph;

/* Graph Constructor */
Graph NewGraph(int size);
EuclidGraph NewEuclid(int size);
MatrixGraph NewMatrix(int size);
Graph CopyGraph(Graph g);

/* Edge Constructor */
void AddEdge(Graph g, int n, int m, int label);
Edge FindEdge(Graph graph, int i, int j);
int RemoveEdge(Graph graph, Edge edge);
int NumEdges(Graph g);

/* Edmond matching  */
int *Weighted_Match(Graph gptr, int type, int maximize);
void Initialize(int maximize);
void FreeUp();

/* Util  */
void POINTER(int u, int v, int e);
void SCAN(int x, int del);
void UNPAIR(int oldbase, int oldmate);
void REMATCH(int firstmate, int e);
void UNLINK(int oldbase);
void SET_BOUNDS();
void UNPAIR_ALL();
void SetUp(Graph gptr, int type);
void SetStandard(Graph graph);
void SetEuclid(EuclidGraph graph);
void SetMatrix(MatrixGraph graph);
void PAIR(int *outcome);
void MERGE_PAIRS(int v);
void LINK_PATH(int e);
void INSERT_PAIR();

/* Euclidean distance routines */
int eucdist(EuclidGraph graph, int i, int j);
int eucdist2(EuclidGraph graph, int i, int j);
int eucdistsq(EuclidGraph graph, int i, int j);


/* Graph I/O routines */
//Graph ReadGraph (int *size, char file[]);
//void WriteGraph (Graph graph,char file[]);
//EuclidGraph ReadEuclid(int *size,char file[]);
//void WriteEuclid(EuclidGraph graph, char file[]);
//MatrixGraph ReadMatrix(int *size, char file[]);
//void WriteMatrix(MatrixGraph graph, char file[]);



#endif /* defined(__Ximulator__edmond__) */
