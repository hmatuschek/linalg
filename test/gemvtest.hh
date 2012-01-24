#ifndef GEMVTEST_HH
#define GEMVTEST_HH

#include "unittest.hh"


class GEMVTest : public UnitTest::TestCase
{
public:
  void testRowMajor();
  void testColMajor();
  void testTransposedRowMajor();
  void testTransposedColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // GEMVTEST_HH
