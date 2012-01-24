#include "gemmtest.hh"
#include "matrix.hh"
#include "blas/gemm.hh"

using namespace Linalg;


void
GEMMTest::testRectRowMajor()
{
  Matrix<double> A(2,3), B(3,2), C(2,2);
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(0,2) = 6;
  B(0,0) = 7; B(0,1) = 8;
  B(1,0) = 9; B(1,1) = 10;
  B(2,0) = 11; B(2,1) = 12;
  C(0,0) = 0; C(0,1) = 0;
  C(1,0) = 0; C(1,1) = 0;

  Blas::gemm(1., A, B, 0., C);

  UT_ASSERT_EQUAL(C(0,0), 58.);
  UT_ASSERT_EQUAL(C(0,1), 64.);
  UT_ASSERT_EQUAL(C(1,0), 139.);
  UT_ASSERT_EQUAL(C(1,1), 154.);
}



void
GEMMTest::testRectColMajor()
{
  Matrix<double> A(2,3, false), B(3,2, false), C(2,2, false);
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(0,2) = 6;
  B(0,0) = 7; B(0,1) = 8;
  B(1,0) = 9; B(1,1) = 10;
  B(2,0) = 11; B(2,1) = 12;
  C(0,0) = 0; C(0,1) = 0;
  C(1,0) = 0; C(1,1) = 0;

  Blas::gemm(1., A, B, 0., C);

  UT_ASSERT_EQUAL(C(0,0), 58.);
  UT_ASSERT_EQUAL(C(0,1), 64.);
  UT_ASSERT_EQUAL(C(1,0), 139.);
  UT_ASSERT_EQUAL(C(1,1), 154.);
}



UnitTest::TestSuite *
GEMMTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::gemm()");

  s->addTest(new UnitTest::TestCaller<GEMMTest>(
               "Blas::gemm(double[m,n], double[n,m], double[m,m]) (row-major)",
               &GEMMTest::testRectRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMMTest>(
               "Blas::gemm(double[m,n], double[n,m], double[m,m]) (col-major)",
               &GEMMTest::testRectColMajor));

  return s;
}
