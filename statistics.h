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

#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <iostream>

namespace oct {

/// /f$\sqrt{\frac{1}{N} \left(\sum_{i=1}^N x_i^2\right) - \overline{x}^2}/f$
/// is used to calculate the standard deviation.
///
/// No checks are made for overflow!  This class should not be used for critical
/// calculations.
template <typename T>
class statistics
{
public:
  statistics() : _N(0), _sum(0), _squares(0) {}
  ~statistics() {}

  void operator()(T sample) {
    add(sample);
  }

  void add(T sample)
  {
    if (_N == 0) {
      _min = sample;
      _max = sample;
    }
    ++_N;
    _sum += sample;
    _squares += (sample * sample);
    _min = (_min<sample)?_min:sample;
    _max = (_max>sample)?_max:sample;
  }

  size_t N() const
  { return _N; }

  T mean() const
  { return _sum / (T) _N; }

  T variance() const
  {
    T avg = mean();
    return _squares / (T) _N - avg * avg;
  }

  T std_dev() const
  { return sqrt(variance()); }

  T min() const
  { return _min; }

  T maximum() const
  { return _max; }

  T range() const {
    return _max - _min;
  }

private:
  size_t _N;
  T _sum;
  T _squares;
  T _min;
  T _max;
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const statistics<T>& stats)
{
  out << "Sample size = " << stats.N() << " Average = " << stats.mean() << " Standard deviation = " << stats.std_dev()
    << " Min = " << stats.min() << " Max = " << stats.maximum();
  return out;
}

}

#endif
