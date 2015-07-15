// test_Range.cxx

// David Adams
// May 2015
//
// Test script for Range.

#include "LbTupler/Range.h"

#include <string>
#include <iostream>
#include <cassert>

using std::string;
using std::cout;
using std::endl;

typedef Range<int> IntRange;
int main() {
  const string myname = "test_Range: ";
  cout << myname << "Starting test" << endl;
  IntRange rng(1,10);
  assert( *rng.begin() == 1 );
  assert( *rng.end() == 11 );
  assert( rng.first() == 1 );
  assert( rng.last() == 10 );
  int count = 0;
  for ( auto i : rng ) {
    cout << myname << "  " << i << endl;
    ++count;
  }
  assert( count == 10 );
  cout << myname << "Ending test" << endl;
  return 0;
}
