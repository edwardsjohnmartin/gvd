#ifndef __OCT_NODE_H__
#define __OCT_NODE_H__

namespace Karras {

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
  void set_child(const int quadrant, const int child) {
    children[quadrant] = child;
    if (child > -1) {
      leaf &= ~leaf_masks[quadrant];
    } else {
      leaf |= leaf_masks[quadrant];
    }
  }
  const int& operator[](const int i) const {
    if (!children) throw std::logic_error("leaf can't get children");
    return children[i];
  }

 private:
  int* children;
  unsigned char leaf;
};

} // namespace

#endif
