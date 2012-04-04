#ifndef ARRAYTEST_HH
#define ARRAYTEST_HH

#include "unittest.hh"

class ArrayTest : public UnitTest::TestCase
{
public:
  void testAssignment();
  void testValueAssignment();

public:
  static UnitTest::TestSuite *suite();
};

#endif // ARRAYTEST_HH
