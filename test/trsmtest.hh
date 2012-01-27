#ifndef TRSMTEST_HH
#define TRSMTEST_HH

#include "unittest.hh"


class TRSMTest : public UnitTest::TestCase
{
public:
  void testUpperRowMajor();
  void testUpperColMajor();


public:
  static UnitTest::TestSuite *suite();
};

#endif // TRSMTEST_HH
