#ifndef TRTRITEST_HH
#define TRTRITEST_HH

#include "unittest.hh"


class TRTRITest : public UnitTest::TestCase
{
public:
  void testUpperRowMajor();
  void testUpperTransRowMajor();
  void testTransLowerRowMajor();

  void testUpperColMajor();
  void testUpperTransColMajor();
  void testTransLowerColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // TRTRITEST_HH
