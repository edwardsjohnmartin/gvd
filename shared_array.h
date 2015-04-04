/*******************************************************
 ** Generalized Voronoi Diagram Project               **
 ** Copyright (c) 2015 John Martin Edwards            **
 ** Scientific Computing and Imaging Institute        **
 ** 72 S Central Campus Drive, Room 3750              **
 ** Salt Lake City, UT 84112                          **
 **                                                   **
 ** For information about this project contact        **
 ** John Edwards at                                   **
 **    edwardsjohnmartin@gmail.com                    **
 ** or visit                                          **
 **    sci.utah.edu/~jedwards/research/gvd/index.html **
 *******************************************************/

#ifndef __SHARED_ARRAY_H__
#define __SHARED_ARRAY_H__

#include <assert.h>

namespace oct {

template <typename T>
class shared_array {
 public:
  shared_array() : _value(0), _refs(0) {}
  explicit shared_array(T* value) : _value(value), _refs(new int(1)) {}
  shared_array(shared_array<T> const & a)
      : _value(a._value), _refs(a._refs) {
    if (_refs) {
      ++(*_refs);
    }
  }
  ~shared_array() {
    dec();
  }
  shared_array& operator=(shared_array const& a) {
    dec();
    _value = a._value;
    _refs = a._refs;
    if (_refs) {
      ++(*_refs);
    }
    return *this;
  }
  T& operator*() { return *_value; }
  const T& operator*() const { return *_value; }
  T* operator->() const {
    return _value;
  }
  operator bool() const { return _value != 0; }
  T& operator[](int i) { return _value[i]; }
  const T& operator[](int i) const { return _value[i]; }

  T* get() const {
    return _value;
  }

  void reset(T* value) {
    dec();
    _value = value;
    _refs = new int(1);
  }

 private:
  void dec() {
    if (_refs) {
      --(*_refs);
      assert(*_refs >= 0);
      if (*_refs == 0) {
        delete[] _value;
        delete _refs;
      }
    }
  }

 private:
  T* _value;
  int* _refs;
};
}

#endif
