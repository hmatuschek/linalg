#ifndef MATRIXTEST_HH
#define MATRIXTEST_HH

#include "unittest.hh"


class MatrixTest : public UnitTest::TestCase
{
public:
  void testRectRowToColMajor();
  void testSquareRowToColMajor();
  void testSubRowToColMajor();

  void testRectColToRowMajor();
  void testSquareColToRowMajor();
  void testSubColToRowMajor();

  void testFromData();

public:
  static UnitTest::TestSuite *suite();
};

#endif // MATRIXTEST_HH
