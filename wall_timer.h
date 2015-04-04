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

#ifndef __WALL_TIMER_H__
#define __WALL_TIMER_H__

#include <time.h>
#include <sstream>
#include <string>
#include <iostream>

namespace oct {

class WallTimer {
 public:
  WallTimer(const std::string& msg, const std::string& id = "",
        bool initial_msg = false)
      // : _msg(msg), _id(id), _t(clock()), _active(true),
      : _msg(msg), _id(id), _active(true),
        _initial_msg(initial_msg), _output(true), _alive(true) {

    gettimeofday(&_t, 0);
    initial();
  }

  ~WallTimer() {
    stop();
  }

  void kill() {
    _alive = false;
  }

  void suspend() {
    throw std::logic_error("not implemented");
    if (!_alive) return;
    _active = false;
    // _t = clock()-_t;
  }

  void resume() {
    throw std::logic_error("not implemented");
    if (!_alive) return;
    _active = true;
    // _t = clock()-_t;
  }

  void restart(const std::string& msg) {
    if (!_alive) return;
    stop();
    _msg = msg;
    // _t = clock();
    gettimeofday(&_t, 0);
    _active = true;
    initial();
  }

  void reset(const std::string& msg) {
    if (!_alive) return;
    restart(msg);
  }

  void stop() {
    if (!_alive) return;
    if (_active) {
      // const double t = (clock()-_t)/static_cast<double>(CLOCKS_PER_SEC);
      timeval end_t;
      gettimeofday(&end_t, 0);
      const double t = (end_t.tv_sec-_t.tv_sec) +
          (end_t.tv_usec-_t.tv_usec)/(double)1000000;
      std::stringstream ss;
      if (_output) {
        if (!_id.empty()) std::cout << _id << " ";
        std::cout << _msg << " " << t << " sec" << std::endl;
      }
      _active = false;
    }
  }

  void set_output(bool output) {
    _output = output;
  }

 private:
  void initial() const {
    if (!_alive) return;
    if (_initial_msg) {
      if (_output) {
        if (!_id.empty()) std::cout << _id << " ";
        std::cout << _msg << "..." << std::endl;
      }
    }
  }

 private:
  std::string _msg;
  std::string _id;
  // time_t _t;
  timeval _t;
  bool _active;
  bool _initial_msg;
  bool _output;
  bool _alive;
};

}  // end namespace

#endif
