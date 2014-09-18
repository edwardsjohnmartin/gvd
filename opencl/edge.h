#ifndef __EDGE_H__
#define __EDGE_H__

typedef struct {
  int s[2];
} Edge;

inline Edge make_edge(int i, int j) {
  int v[2] = { i, j };
  return *(Edge*)v;
}

#endif
