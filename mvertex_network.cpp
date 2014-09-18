#include "./mvertex_network.h"

NAMESPACE_OCT_BEGIN

MVertexNetwork make_mvertex_network() {
  MVertexNetwork mvn;
  mvn.header.reset(new int[2]);
  SetNumVertices(0, &mvn);
  mvn.vertex_array_capacity = 64;
  mvn.vertices = shared_array<Vertex>(new Vertex[64]);

  SetNumCPoints(0, &mvn);
  mvn.cpoint_array_capacity = 0;
  mvn.cpoints = shared_array<GeomPoint>();

  UVertexNetwork uvn = make_vertex_network(mvn);
  Initialize(&uvn);
  update_mvertex_network(uvn, mvn);
  return mvn;
}

UVertexNetwork make_vertex_network(MVertexNetwork mvn) {
  UVertexNetwork vn;
  vn.header = mvn.header.get();
  // vn.num_vertices = NumVertices(mvn);
  vn.vertices = mvn.vertices.get();

  // vn.num_cpoints = NumCPoints(mvn);
  vn.cpoints = mvn.cpoints.get();

  return vn;
}

void update_mvertex_network(UVertexNetwork vn, MVertexNetwork& mvn) {
  assert(vn.vertices == mvn.vertices.get());
  assert(vn.cpoints == mvn.cpoints.get());

  // mvn.num_vertices = NumVertices(vn);
  // mvn.num_cpoints = NumCPoints(vn);
  SetNumVertices(NumVertices(vn), &mvn);
  SetNumCPoints(NumCPoints(vn), &mvn);
}

NAMESPACE_OCT_END

