#include "trmmtest.hh"
#include "trimatrix.hh"
#include "blas/trmm.hh"

using namespace Linalg;


void
TRMMTest::testUpperRowMajor()
{
  Matrix<double> A(2,2);
  Matrix<double> B(2,3);

  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;
  B(0,0) = 5; B(0,1) = 6; B(0,2) = 7;
  B(1,0) = 8; B(1,1) = 9; B(1,2) = 10;

  // Multiply with upper triangular of A
  Matrix<double> tmp(B.copy());
  Blas::trmm(true, 1, triu(A), tmp);
  UT_ASSERT_EQUAL(tmp(0,0), 21.);  UT_ASSERT_EQUAL(tmp(0,1), 24.);  UT_ASSERT_EQUAL(tmp(0,2), 27.);
  UT_ASSERT_EQUAL(tmp(1,0), 32.);  UT_ASSERT_EQUAL(tmp(1,1), 36.);  UT_ASSERT_EQUAL(tmp(1,2), 40.);

  tmp = B.copy();
  Blas::trmm(true, 1, triu(A, true), tmp);
  UT_ASSERT_EQUAL(tmp(0,0), 21.);  UT_ASSERT_EQUAL(tmp(0,1), 24.);  UT_ASSERT_EQUAL(tmp(0,2), 27.);
  UT_ASSERT_EQUAL(tmp(1,0), 8.);   UT_ASSERT_EQUAL(tmp(1,1), 9.);   UT_ASSERT_EQUAL(tmp(1,2), 10.);
}


void
TRMMTest::testUpperTransRowMajor()
{
  Matrix<double> A(2,2);
  Matrix<double> B(2,3);

  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;
  B(0,0) = 5; B(0,1) = 6; B(0,2) = 7;
  B(1,0) = 8; B(1,1) = 9; B(1,2) = 10;

  // Multiply with upper triangular of A
  Matrix<double> tmp(B.copy());
  Blas::trmm(true, 1, triu(A).t(), tmp);
  UT_ASSERT_EQUAL(tmp(0,0),  5.);  UT_ASSERT_EQUAL(tmp(0,1),  6.);  UT_ASSERT_EQUAL(tmp(0,2),  7.);
  UT_ASSERT_EQUAL(tmp(1,0), 42.);  UT_ASSERT_EQUAL(tmp(1,1), 48.);  UT_ASSERT_EQUAL(tmp(1,2), 54.);

  tmp = B.copy();
  Blas::trmm(true, 1, triu(A, true).t(), tmp);
  UT_ASSERT_EQUAL(tmp(0,0),  5.);  UT_ASSERT_EQUAL(tmp(0,1),  6.);  UT_ASSERT_EQUAL(tmp(0,2),  7.);
  UT_ASSERT_EQUAL(tmp(1,0), 18.);  UT_ASSERT_EQUAL(tmp(1,1), 21.);  UT_ASSERT_EQUAL(tmp(1,2), 24.);
}


void
TRMMTest::testTransUpperRowMajor()
{
  Matrix<double> A(2,2);
  Matrix<double> B(2,3);

  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;
  B(0,0) = 5; B(0,1) = 6; B(0,2) = 7;
  B(1,0) = 8; B(1,1) = 9; B(1,2) = 10;

  // Multiply with upper triangular of A
  Matrix<double> tmp(B.copy());
  Blas::trmm(true, 1, triu(A.t()), tmp);
  UT_ASSERT_EQUAL(tmp(0,0), 29.);  UT_ASSERT_EQUAL(tmp(0,1), 33.);  UT_ASSERT_EQUAL(tmp(0,2), 37.);
  UT_ASSERT_EQUAL(tmp(1,0), 32.);  UT_ASSERT_EQUAL(tmp(1,1), 36.);  UT_ASSERT_EQUAL(tmp(1,2), 40.);

  tmp = B.copy();
  Blas::trmm(true, 1, triu(A.t(), true), tmp);
  UT_ASSERT_EQUAL(tmp(0,0), 29.);  UT_ASSERT_EQUAL(tmp(0,1), 33.);  UT_ASSERT_EQUAL(tmp(0,2), 37.);
  UT_ASSERT_EQUAL(tmp(1,0), 8.);   UT_ASSERT_EQUAL(tmp(1,1), 9.);   UT_ASSERT_EQUAL(tmp(1,2), 10.);
}



void
TRMMTest::testUpperColMajor()
{
  Matrix<double> A(2,2, false);
  Matrix<double> B(2,3, false);

  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;
  B(0,0) = 5; B(0,1) = 6; B(0,2) = 7;
  B(1,0) = 8; B(1,1) = 9; B(1,2) = 10;

  // Multiply with upper triangular of A
  Matrix<double> tmp(B.copy());
  Blas::trmm(true, 1, triu(A), tmp);
  UT_ASSERT_EQUAL(tmp(0,0), 21.);  UT_ASSERT_EQUAL(tmp(0,1), 24.);  UT_ASSERT_EQUAL(tmp(0,2), 27.);
  UT_ASSERT_EQUAL(tmp(1,0), 32.);  UT_ASSERT_EQUAL(tmp(1,1), 36.);  UT_ASSERT_EQUAL(tmp(1,2), 40.);

  // Multiply with upper triangular of A with unit diagonal
  tmp = B.copy();
  Blas::trmm(true, 1, triu(A, true), tmp);
  UT_ASSERT_EQUAL(tmp(0,0), 21.);  UT_ASSERT_EQUAL(tmp(0,1), 24.);  UT_ASSERT_EQUAL(tmp(0,2), 27.);
  UT_ASSERT_EQUAL(tmp(1,0), 8.);   UT_ASSERT_EQUAL(tmp(1,1), 9.);   UT_ASSERT_EQUAL(tmp(1,2), 10.);
}



UnitTest::TestSuite *
TRMMTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::trmm()");

  s->addTest(new UnitTest::TestCaller<TRMMTest>(
               "Blas::trmm(left, triu(double[m,m]), double[m,n]) (row-major)",
               &TRMMTest::testUpperRowMajor));

  s->addTest(new UnitTest::TestCaller<TRMMTest>(
               "Blas::trmm(left, triu(double[m,m])::t(), double[m,n]) (row-major)",
               &TRMMTest::testUpperTransRowMajor));

  s->addTest(new UnitTest::TestCaller<TRMMTest>(
               "Blas::trmm(left, triu(double[m,m]::t()), double[m,n]) (row-major)",
               &TRMMTest::testTransUpperRowMajor));

  s->addTest(new UnitTest::TestCaller<TRMMTest>(
               "Blas::trmm(left, triu(double[m,m]), double[m,n]) (col-major)",
               &TRMMTest::testUpperColMajor));

  return s;
}
