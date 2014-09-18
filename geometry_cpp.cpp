//------------------------------------------------------------
// Test code to be used only in c++
//------------------------------------------------------------

#include "./geometry_cpp.h"
#include "./opencl/tribox.h"

#include "./vector3.h"
#include "./triangle_cpp.h"

using namespace std;

NAMESPACE_OCT_BEGIN

LabeledGeometry3 Convert(Geometry geometry) {
  using namespace std;
  // const int label = geometry.label;
  const int label = g_label(geometry);
  // const int m = geometry.n;
  const int m = g_m(geometry);
  vector<Face> faces;
  for (int k = 0; k < m; ++k) {
    faces.push_back(geometry.faces[k]);
  }
  LabeledGeometry3 lg(faces, label);
  return lg;
}

std::vector<LabeledGeometry3> Convert(Geometries geometries) {
  using namespace std;
  vector<LabeledGeometry3> vcell_geoms;
  // const int n = geometries.n;
  const int n = g_n(geometries);
  for (int j = 0; j < n; ++j) {
    const Geometry geometry = get_geometry(j, geometries);
    LabeledGeometry3 lg = Convert(geometry);
    if (!lg.GetTriangles().empty()) {
      vcell_geoms.push_back(lg);
    }
  }
  return vcell_geoms;
}

std::vector<std::vector<LabeledGeometry3> > Convert(
    Vi2Geometries vi2geometries) {
  using namespace std;
  // Copy arrays to base2goemetries2
  vector<vector<LabeledGeometry3> > base2geometries;
  for (int vi = 0; vi < g_N(vi2geometries); ++vi) {
    // const Geometries geometries =
    //     get_geometries(vi, vi2geometries);
    vector<LabeledGeometry3> vcell_geoms;
    if (!gcell_is_empty(vi, vi2geometries)) {
      const Geometries geometries =
          get_geometries(vi, vi2geometries);
      vcell_geoms = Convert(geometries);
    }
    base2geometries.push_back(vcell_geoms);
  }
  return base2geometries;
}

// The size of an array containing geometries belonging to a
// single octree cell.
size_t lgsize(const LabeledGeometry3& geom) {
  if (geom.size() == 0) return 0;
  // 2: label, size
  return 2 + geom.size() * DIM;
}

size_t lgsize(const std::vector<LabeledGeometry3>& geoms) {
  int gsize = 0;
  for (int i = 0; i < geoms.size(); ++i) {
    gsize += lgsize(geoms[i]);
  }
  if (gsize == 0) return 0;
  return geoms.size() + gsize;
}

size_t lgsize(const std::vector<std::vector<LabeledGeometry3> >& geoms) {
  int geom_size_sum = 0;
  for (int i = 0; i < geoms.size(); ++i) {
    geom_size_sum += lgsize(geoms[i]);
  }
  return geoms.size() + geom_size_sum;
}

// returns the number of bytes written
int Convert(
    int* array,
    const LabeledGeometry3& g) {
  int offset = 0;
  array[offset++] = g.GetLabel();
  array[offset++] = g.size();
  for (int i = 0; i < g.size(); ++i) {
    const Face& face = g.GetPrimitives()[i];
    for (int k = 0; k < DIM; ++k) {
      array[offset++] = face.s[k];
    }
  }
  return offset;
}

// returns the number of bytes written
int Convert(
    int* array,
    const std::vector<LabeledGeometry3>& g3) {
  int offset = g3.size();
  for (int i = 0; i < g3.size(); ++i) {
    array[i] = offset;
    offset += Convert(array+offset, g3[i]);
  }
  return offset;
}

// shared_array<int> Convert(
//     const std::vector<std::vector<LabeledGeometry3> >& lg,
//     Vi2Geometries& vi2geometries) {
//   const int size = lgsize(lg);
//   shared_array<int> array(new int[size]);
//   int offset = lg.size();
//   for (int vi = 0; vi < lg.size(); ++vi) {
//     const vector<LabeledGeometry3>& geoms = lg[vi];
//     if (geoms.empty()) {
//       array[vi] = -1;
//     } else {
//       array[vi] = offset;
//       offset += Convert(array.get()+offset, geoms);
//     }
//   }
//   vi2geometries = make_vi2geometries(array.get(), size);
//   return array;
// }

// Write geometries to vi2geometries
shared_array<int> Convert(
    const int num_vertices,
    const vector<LabeledGeometry3>& geometries,
    Vi2Geometries& vi2geometries) {

  // Count array size
  int size = vi2g_offset + num_vertices + lgsize(geometries);

  // Populate the array
  shared_array<int> array(new int[size]);
  array[0] = size;
  array[1] = size;
  // array[0] = num_vertices;
  array[vi2g_offset] = vi2g_offset + num_vertices;
  for (int vi = 1; vi < num_vertices; ++vi) {
    array[vi2g_offset+vi] = -1;
  }
  Convert(array.get() + vi2g_offset + num_vertices, geometries);
  // vi2geometries = make_vi2geometries(array.get(), size);
  vi2geometries = make_vi2geometries(array.get());

  return array;
}

// // The array size in geometries must be large enough to handle old_geoms.
// void Insert(const std::vector<LabeledGeometry3>& old_geometries,
//                    const int vi, Vi2Geometries& vi2geometries) {
//   using namespace std;

//   bool update_free_offset = false;
//   if (gcell_is_empty(vi, vi2geometries)) {
//     update_free_offset = true;
//     vi2geometries.offsets[vi] = vi2geometries.free_offset;
//   }

//   Geometries geometries = get_geometries(vi, vi2geometries);

//   const int n = old_geometries.size();
//   g_set_n(n, geometries);
//   int offset = n;
//   for (int i = 0; i < n; ++i) {
//     geometries.offsets[i] = offset;
//     const LabeledGeometry3& old_geometry = old_geometries[i];
//     Geometry geometry = get_geometry(i, geometries);
//     vector<Face> faces = old_geometry.GetPrimitives();
//     const int label = old_geometry.GetLabel();
//     const int m = faces.size();
//     assert(m);
//     // geometry.label = label;
//     g_set_label(label, geometry);
//     // geometry.n = m;
//     g_set_m(m, geometry);
//     for (int j = 0; j < m; ++j) {
//       geometry.faces[j] = faces[j];
//     }
//     // increment offset
//     offset += geometry_size(geometry);
//   }

//   if (update_free_offset) {
//     vi2geometries.free_offset += geometries_size(geometries);
//   }
// }

std::ostream& operator<<(
    std::ostream& out, Vi2Geometries vi2geometries) {
  using namespace std;
  const int N = g_N(vi2geometries);
  out << " num vertices = " << N
      // << " array size = " << vi2geometries.array_size
      // << " free offset = " << vi2geometries.free_offset
      << " array size = " << g_array_size(vi2geometries)
      << " free offset = " << g_free_offset(vi2geometries)
      << endl;
  for (int vi = 0; vi < N; ++vi) {
    Geometries geometries = get_geometries(vi, vi2geometries);
    const int n = g_n(geometries);
    if (n > 0) {
      // See if all geometries are empty
      bool empty = true;
      for (int j = 0; empty && j < n; ++j) {
        const Geometry geometry = get_geometry(j, geometries);
        const int m = g_m(geometry);
        empty = (m == 0);
      }

      if (empty)
        out << "* ";
      else
        out << "  ";
      out << "vi = " << vi
          // << " array offset = " << vi2geometries.offsets[vi]
          << " array = " << vi2geometries.offsets[vi]
          << "-" << vi2geometries.offsets[vi]+geometries_size(geometries)
          << " asize = " << geometries_size(geometries)
          << " num_geoms = " << n << " (";
      // geometry labels
      for (int gi = 0; gi < n; ++gi) {
        out << g_label(get_geometry(gi, geometries));
        if (gi < n-1)
          out << " ";
      }
      out << ") (";
      // geometry triangle counts
      for (int gi = 0; gi < n; ++gi) {
        out << g_m(get_geometry(gi, geometries));
        if (gi < n-1)
          out << " ";
      }
      out << ") (";
      // "checksum" of triangle indices
      for (int gi = 0; gi < n; ++gi) {
        Geometry geometry = get_geometry(gi, geometries);
        const int m = g_m(geometry);
        int sum = 0;
        for (int ti = 0; ti < m; ++ti) {
          Triangle t = geometry.faces[ti];
          for (int k = 0; k < 3; ++k) {
            sum += t.s[k];
          }
        }
        out << sum;
        if (gi < n-1)
          out << " ";
      }
      out << ")" << endl;
    }
  }
  return out;
}

void Compare(
    const std::vector<std::vector<LabeledGeometry3> >& base2geometries,
    Vi2Geometries vi2geometries) {
  using namespace std;
  for (int vi = 0; vi < base2geometries.size(); ++vi) {
    Geometries geometries = get_geometries(vi, vi2geometries);
    CompareAndExit((int)base2geometries[vi].size(), g_n(geometries),
                   vi, "number of geometries");

    const int n = g_n(geometries);
    for (int j = 0; j < n; ++j) {
      const Geometry geometry = get_geometry(j, geometries);
      CompareAndExit((int)base2geometries[vi][j].size(), g_m(geometry),
                     vi, "number of primitives");
      CompareAndExit(base2geometries[vi][j].GetLabel(), g_label(geometry),
                     vi, "label");
      for (int k = 0; k < g_m(geometry); ++k) {
        Face f0 = base2geometries[vi][j].GetPrimitives()[k];
        Face f1 = geometry.faces[k];
        CompareAndExit(f0, f1, vi, "faces");
      }
    }
  }
}

NAMESPACE_OCT_END
