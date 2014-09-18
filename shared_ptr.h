#ifndef __SHARED_PTR_H__
#define __SHARED_PTR_H__

#include <assert.h>

namespace oct {

template <typename T>
class shared_ptr {
 public:
  shared_ptr() : _value(0), _refs(0) {}
  explicit shared_ptr(T* value) : _value(value), _refs(new int(1)) {}
  shared_ptr(shared_ptr<T> const & a)
      : _value(a._value), _refs(a._refs) {
    if (_refs) {
      ++(*_refs);
    }
  }
  ~shared_ptr() {
    dec();
  }
  shared_ptr& operator=(shared_ptr const& a) {
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
        delete _value;
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
