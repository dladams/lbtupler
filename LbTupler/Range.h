// Range.h

#ifndef Range_H
#define Range_H

// David Adams
// May 2015
//
// Class to hold and iterate over a range of values.
// Can be used in range-based for loops, e.g. to count from 1 to 10:
//   Range<int> counter(1,10);
//   for ( int count : counter ) {...}

template<class T>
class Range {
public:
  Range(const T& val) : m_val1(val), m_val2(val) { };
  Range(const T& val1, const T& val2) : m_val1(val1), m_val2(val2) { };
  T& first() { return m_val1; }
  T& last() { return m_val2; }
  const T& first() const { return m_val1; }
  const T& last() const { return m_val2; }
  const T& begin() const { return m_val1; }
  const T& end() const { T val = m_val2; return ++val; }
private:
  T m_val1;
  T m_val2;
};

#endif
