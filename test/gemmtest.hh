#ifndef GEMMTEST_HH
#define GEMMTEST_HH

#include "unittest.hh"


class GEMMTest : public UnitTest::TestCase
{
public:
  void testRectRowMajor();
  void testRectColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // GEMMTEST_HH
