#ifndef TRIMATRIXTEST_HH
#define TRIMATRIXTEST_HH

#include "unittest.hh"


class TriMatrixTest : public UnitTest::TestCase
{
public:
  void testRowToColMajor();

public:
  UnitTest::TestSuite *suite();
};

#endif // TRIMATRIXTEST_HH
