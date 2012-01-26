#ifndef TRMVTEST_HH
#define TRMVTEST_HH

#include "unittest.hh"
#include "matrix.hh"
#include "vector.hh"


class TRMVTest : public UnitTest::TestCase
{
private:
  Linalg::Matrix<double> Arow;
  Linalg::Matrix<double> Acol;
  Linalg::Vector<double> x;

public:
  virtual void setUp();

public:
  void testUpperRowMajor();
  void testUpperTransRowMajor();
  void testTransUpperRowMajor();

  void testUpperColMajor();
  void testUpperTransColMajor();
  void testTransUpperColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // TRMVTEST_HH
