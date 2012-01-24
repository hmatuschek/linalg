#ifndef DOTTEST_HH
#define DOTTEST_HH

#include "unittest.hh"


class DOTTest : public UnitTest::TestCase
{
public:
  void testVectorVector();

public:
  static UnitTest::TestSuite *suite();
};

#endif // DOTTEST_HH
