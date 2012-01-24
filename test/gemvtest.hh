#ifndef GEMVTEST_HH
#define GEMVTEST_HH

#include "unittest.hh"


class GEMVTest : public UnitTest::TestCase
{
public:
  void testSquareRowMajor();
  void testSquareColMajor();
  void testSquareTransposedRowMajor();
  void testSquareTransposedColMajor();

  void testRectRowMajor();
  void testRectColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // GEMVTEST_HH
