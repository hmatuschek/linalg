#include "trtrstest.hh"
#include "matrix.hh"
#include "trimatrix.hh"
#include "lapack/trtrs.hh"

using namespace Linalg;


void
TRTRSTest::testUpperRowMajor()
{
  Matrix<double> A(3,3);
  Matrix<double> B(3,2);

  A(0,0)=1; A(0,1)=2; A(0,2)=3;
  A(1,0)=4; A(1,1)=5; A(1,2)=6;
  A(2,0)=7; A(2,1)=8; A(2,2)=9;

  B(0,0)=10; B(0,1)=11;
  B(1,0)=12; B(1,1)=13;
  B(2,0)=14; B(2,1)=15;

  TriMatrix<double> Au(triu(A));
  Matrix<double> Bc(B.copy());
  Lapack::trtrs<>(Au, Bc);

  UT_ASSERT_NEAR(Bc(0,0), 42./10 + 2./30);
  UT_ASSERT_NEAR(Bc(0,1), 48./10);
  UT_ASSERT_NEAR(Bc(1,0), 5./10 + 1./30);
  UT_ASSERT_NEAR(Bc(1,1), 6./10);
  UT_ASSERT_NEAR(Bc(2,0), 1. + 5./9);
  UT_ASSERT_NEAR(Bc(2,1), 1. + 2./3);
}



UnitTest::TestSuite *
TRTRSTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Lapack::trtrs()");

  s->addTest(new UnitTest::TestCaller<TRTRSTest>(
               "Lapack::trtrs(triu(double[m,m]), double[m,n]) (row-major)",
               &TRTRSTest::testUpperRowMajor));

  return s;
}
