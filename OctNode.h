#ifndef __OCT_NODE_H__
#define __OCT_NODE_H__

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

 private:
  int* children;
  unsigned char leaf;
};


struct OctCell {
  OctCell() : parent(0) {}
  OctCell(const intn origin_, const int width_,
          OctNode const* parent_, int octant_,
          OctNode const* node_, int data_)
      : origin(origin_), width(width_),
        parent(parent_), octant(octant_),
        node(node_), data(data_) {}

  intn get_origin() const { return origin; }
  int get_width() const { return width; }
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

 private:
  intn origin;
  int width;
  OctNode const* parent;
  int octant;
  OctNode const* node;
  int data;
};


} // namespace

#endif
