#ifndef TRMMTEST_HH
#define TRMMTEST_HH

#include "unittest.hh"


class TRMMTest : public UnitTest::TestCase
{
public:
  void testRowMajor();
  void testTransposedRowMajor();

  void testColMajor();
  void testTransposedColMajor();

public:
  static UnitTest::TestSuite *suite();
};


#endif // TRMMTEST_HH
