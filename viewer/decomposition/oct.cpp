const int DIM = 2; // TODO

VertexNetwork ExtendOctree(
    VertexNetwork tree_verts,
    vector<float2> obj_verts)
{
    ;
}



VertexNetwork BuildOctree(
    const vector<vector<floatn> > &label2gverts,    // list of vertices
    const vector<vector<Face> > &label2faces,       // list of edges
    const BoundingBox<floatn> &bb,                  // root cell box region
    const OctreeOptions &o)                         // options object
{
    // Establish dimension... why not just use DIM?
    static const int D = DIM;

    // Sanity check for correct input lists.
    if(label2faces.size() != label2gverts.size())
        return VertexNetwork();
    if(label2faces.empty())
        return VertexNetwork();
    if(label2faces[0].empty())
        return VertexNetwork();

    const int num_labels = label2gverts.size();

    Timer t("Convert vertices", "*BuildOctree*");
    t.set_output(o.timings || o.report_statistics);

    // Find bounding box of vertices
    BoundingBox<floatn> vert_bb;
    for (int j = 0; j < num_labels; ++j) {
        const std::vector<floatn>& fvertices = label2gverts[j];
        for (int i = 0; i < fvertices.size(); ++i) {
            vert_bb(fvertices[i]);
        }
    }

    const float dwidth = bb.size().s[0];

    if (dwidth == 0) return VertexNetwork();

    // Count total number of geometry vertices
    int num_gverts = 0;
    for (int j = 0; j < num_labels; ++j) {
        num_gverts += label2gverts[j].size();
    }

    // Convert vertices to integer coordinates
    std::vector<LabeledGeometry> lgeometries;
    shared_array<intn> gverts(new intn[num_gverts]);
    shared_array<int> gvert_offsets(new int[num_labels]);
    GeomVertices geom_vertices = {
        num_gverts, gverts.get(), label2gverts.size(), gvert_offsets.get() };
    int geom_vert_idx = 0;
    for (int j = 0; j < label2gverts.size(); ++j) {
        gvert_offsets[j] = geom_vert_idx;
        const std::vector<floatn>& fvertices = label2gverts[j];
        const int n = fvertices.size();
        for (int i = 0; i < n; ++i) {
            intn p = make_intn(0);
            for (int k = 0; k < D; ++k) {
                const double d =
                    (kWidth-1) * ((fvertices[i].s[k] - bb.min().s[k]) / dwidth);
                int v = static_cast<int>(d+0.5);
                if (v < 0) {
                    cerr << "Coordinate in dimension " << k << " is less than zero.  d = "
                        << d << " v = " << v << endl;
                    cerr << "  fvertices[i][k] = " << fvertices[i].s[k]
                        << " bb.min()[k] = " << bb.min().s[k] << endl;
                    cerr << "  dwidth = " << dwidth << " kwidth = " << kWidth << endl;
                    v = 0;
                }
                p.s[k] = v;
            }
            gverts[geom_vert_idx++] = p;
            // cout << "Adding geometry vertex: " << p << endl;
        }
        if (n > 0) {
#ifdef OCT3D
            LabeledGeometry lg(label2faces[j], j);
#else
            int2* lgverts = get_geom_vertices(j, geom_vertices);
            LabeledGeometry lg(lgverts, label2faces[j], j);
#endif
            lgeometries.push_back(lg);
        }
    }

    VertexNetwork vertices;
    if (vertices.size() > (1<<D)) {
        cerr << vertices.size() << endl;
        throw logic_error("vertices passed to BuildOctree must be empty");
    }

    if (o.gpu) {
        GpuTest(lgeometries, geom_vertices, bb, vertices, o);

        MVertexNetwork mvertices = make_mvertex_network();
        BuildOctreeGpu(lgeometries, geom_vertices, mvertices, o);
        vertices = VertexNetwork(
                NumVertices(mvertices), mvertices.vertices.get(),
                NumCPoints(mvertices), mvertices.cpoints.get());
    } else {
        BuildOctreeCpu(lgeometries, geom_vertices, bb, vertices, o);
    }

    // return make_mvertex_network(vertices);
    return vertices;
}
