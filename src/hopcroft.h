#ifndef hopcroftMatching_h
#define hopcroftMatching_h

#define MAXN1  50000
#define MAXN2 50000
#define MAXM 400000

void hopcroftInit(int, int);
void hopcroftAddEdge(int u, int v);
void bfs();
bool dfs(int u);
int hopcroft(int *match_result);
#endif /*!hopcroftMatching_h*/