#ifndef TRTRITEST_HH
#define TRTRITEST_HH

#include "unittest.hh"
#include "matrix.hh"


class TRTRITest : public UnitTest::TestCase
{
protected:
  Linalg::Matrix<double> A, Bu, Bl;

public:
  void setUp();

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
