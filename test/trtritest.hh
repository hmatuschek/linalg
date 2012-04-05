#ifndef TRTRITEST_HH
#define TRTRITEST_HH

#include "unittest.hh"
#include "matrix.hh"


class TRTRITest : public UnitTest::TestCase
{
protected:
  Linalg::Matrix<double> cA, cBu, cBl;
  Linalg::Matrix<double> fA, fBu, fBl;

public:
  void setUp();

  void testUpperRowMajor();
  void testLowerRowMajor();

  void testUpperTransRowMajor();
  void testTransLowerRowMajor();

  void testUpperColMajor();
  void testUpperTransColMajor();
  void testTransLowerColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // TRTRITEST_HH
