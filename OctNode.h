#ifndef __OCT_NODE_H__
#define __OCT_NODE_H__

#include "./resln.h"
#include "./bb.h"

namespace Karras {

// An octree node is an internal node of the octree. An octree cell
// is a general term that refers to both internal nodes and leaves.

static const int leaf_masks[] = { 1, 2, 4, 8 };

struct OctNode {
 public:
  OctNode() : children(new int[1<<DIM]), leaf(15) {
    std::fill(children, children + (1<<DIM), -1);
  }
  ~OctNode() {
    delete [] children;
  }
  bool is_leaf(const int i) const {
    return leaf & leaf_masks[i];
  }
  void set_child(const int octant, const int child) {
    children[octant] = child;
    if (child > -1) {
      leaf &= ~leaf_masks[octant];
    } else {
      leaf |= leaf_masks[octant];
    }
  }
  void set_data(const int octant, const int data) {
    if (!is_leaf(octant))
      throw std::logic_error("Trying to set data on a non-leaf cell");
    children[octant] = data;
  }
  const int& operator[](const int i) const {
    if (!children) throw std::logic_error("leaf can't get children");
    return children[i];
  }

  friend std::ostream& operator<<(std::ostream& out, const OctNode& node);
  friend std::istream& operator>>(std::istream& in, OctNode& node);

 private:
  int* children;
  unsigned char leaf;
};

inline std::ostream& operator<<(std::ostream& out, const OctNode& node) {
  for (int i = 0; i < 1<<DIM; ++i) {
    out << node.children[i] << " ";
  }
  out << static_cast<int>(node.leaf);
  return out;
}

inline std::istream& operator>>(std::istream& in, OctNode& node) {
  for (int i = 0; i < 1<<DIM; ++i) {
    in >> node.children[i];
  }
  int leaf;
  in >> leaf;
  node.leaf = static_cast<unsigned char>(leaf);
  return in;
}

inline std::ostream& operator<<(
    std::ostream& out, const std::vector<OctNode>& octree) {
  out << octree.size() << " ";
  for (const OctNode& node : octree) {
    out << node << " ";
  }
  return out;
}

inline std::istream& operator>>(
    std::istream& in, std::vector<OctNode>& octree) {
  int n;
  in >> n;
  octree.resize(n);
  for (int i = 0; i < n; ++i) {
    in >> octree[i];
  }
  return in;
}

struct OctCell {
  OctCell() : parent(0) {}
  OctCell(const intn origin_, const int width_,
          const int parent_idx_,
          OctNode const* parent_, int octant_,
          OctNode const* node_, int data_)
      : origin(origin_), width(width_),
        parent_idx(parent_idx_), parent(parent_), octant(octant_),
        node(node_), data(data_) {}

  intn get_origin() const { return origin; }
  int get_width() const { return width; }
  int get_parent_idx() const { return parent_idx; }
  OctNode const* get_parent() const { return parent; }
  int get_octant() const { return octant; }
  bool is_leaf() const { return parent->is_leaf(octant); }
  OctNode const* get_node() const {
    if (is_leaf()) {
      throw std::logic_error("Cannot get node from a non-leaf cell");
    }
    return node;
  }
  int get_data() const {
    if (!is_leaf()) {
      throw std::logic_error("Cannot get data from a non-leaf cell");
    }
    return (*parent)[octant];
  }
  int get_level(const Resln& resln) {
    int level = 0;
    while ((width << level) < resln.width)
      ++level;
    return level;
  }
  BoundingBox<intn> bb() const {
    return BoundingBox<intn>(origin, origin+make_uni_intn(width));
  }

  friend std::ostream& operator<<(std::ostream& out, const OctCell& cell);
  friend std::istream& operator>>(std::istream& in, OctCell& cell);

 private:
  intn origin;
  int width;
  int parent_idx;
  OctNode const* parent;
  int octant;
  OctNode const* node;
  int data;
};

inline std::ostream& operator<<(std::ostream& out, const OctCell& cell) {
  out << cell.origin << " " << cell.width << " " << cell.parent_idx
      << " " << cell.octant << " " << cell.data;
  return out;
}

inline std::istream& operator>>(std::istream& in, OctCell& cell) {
  in >> cell.origin >> cell.width >> cell.parent_idx
     >> cell.octant >> cell.data;
  return in;
}


} // namespace

#endif
