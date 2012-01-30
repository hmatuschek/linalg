#ifndef TRTRSTEST_HH
#define TRTRSTEST_HH

#include "unittest.hh"


class TRTRSTest : public UnitTest::TestCase
{
public:
  void testUpperRowMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // TRTRSTEST_HH
