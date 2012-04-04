#include "gemmtest.hh"
#include "matrix.hh"
#include "blas/gemm.hh"

using namespace Linalg;


void
GEMMTest::testRectRowMajor()
{
  size_t M = 2; size_t N = 3; size_t K=4;
  Matrix<double> A(M,N), B(N,K), C(M,K);

  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;

  B(0,0) = 1; B(0,1) = 2; B(0,2) = 3; B(0,3) = 4;
  B(1,0) = 5; B(1,1) = 6; B(1,2) = 7; B(1,3) = 8;
  B(2,0) = 9; B(2,1) =10; B(2,2) =11; B(2,3) =12;

  C(0,0) = 1; C(0,1) = 2; C(0,2) = 3; C(0,3) = 4;
  C(1,0) = 5; C(1,1) = 6; C(1,2) = 7; C(1,3) = 8;

  Blas::gemm(1., A, B, 1., C);

  UT_ASSERT_EQUAL(C(0,0), 39.); UT_ASSERT_EQUAL(C(0,1), 46.); UT_ASSERT_EQUAL(C(0,2), 53.); UT_ASSERT_EQUAL(C(0,3), 60.);
  UT_ASSERT_EQUAL(C(1,0), 88.); UT_ASSERT_EQUAL(C(1,1),104.); UT_ASSERT_EQUAL(C(1,2),120.); UT_ASSERT_EQUAL(C(1,3),136.);
}


void
GEMMTest::testRectColMajor()
{
  size_t M = 2; size_t N = 3; size_t K=4;
  Matrix<double> A(M,N,false), B(N,K,false), C(M,K,false);

  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;

  B(0,0) = 1; B(0,1) = 2; B(0,2) = 3; B(0,3) = 4;
  B(1,0) = 5; B(1,1) = 6; B(1,2) = 7; B(1,3) = 8;
  B(2,0) = 9; B(2,1) =10; B(2,2) =11; B(2,3) =12;

  C(0,0) = 1; C(0,1) = 2; C(0,2) = 3; C(0,3) = 4;
  C(1,0) = 5; C(1,1) = 6; C(1,2) = 7; C(1,3) = 8;

  Blas::gemm(1., A, B, 1., C);

  UT_ASSERT_EQUAL(C(0,0), 39.); UT_ASSERT_EQUAL(C(0,1), 46.); UT_ASSERT_EQUAL(C(0,2), 53.); UT_ASSERT_EQUAL(C(0,3), 60.);
  UT_ASSERT_EQUAL(C(1,0), 88.); UT_ASSERT_EQUAL(C(1,1),104.); UT_ASSERT_EQUAL(C(1,2),120.); UT_ASSERT_EQUAL(C(1,3),136.);
}


void
GEMMTest::testRectTransposedRowMajor()
{
  Matrix<double> A(3,2), B(3,2), C(2,2);

  A(0,0) = 1;  A(0,1) = 2;
  A(1,0) = 3;  A(1,1) = 4;
  A(2,0) = 5;  A(2,1) = 6;

  B(0,0) = 7;  B(0,1) = 8;
  B(1,0) = 9;  B(1,1) = 10;
  B(2,0) = 11; B(2,1) = 12;

  C(0,0) = 0;  C(0,1) = 0;
  C(1,0) = 0;  C(1,1) = 0;

  Blas::gemm(1., A.t(), B, 0., C);

  UT_ASSERT_EQUAL(C(0,0), 89.);
  UT_ASSERT_EQUAL(C(0,1), 98.);
  UT_ASSERT_EQUAL(C(1,0), 116.);
  UT_ASSERT_EQUAL(C(1,1), 128.);
}


void
GEMMTest::testSquareRowMajor()
{
  Matrix<double> A(2,2), B(2,2), C(2,2);
  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;
  B(0,0) = 5; B(0,1) = 6;
  B(1,0) = 7; B(1,1) = 8;
  C(0,0) = 0; C(0,1) = 0;
  C(1,0) = 0; C(1,1) = 0;

  Blas::gemm(1., A, B, 0., C);

  UT_ASSERT_EQUAL(C(0,0), 19.);
  UT_ASSERT_EQUAL(C(0,1), 22.);
  UT_ASSERT_EQUAL(C(1,0), 43.);
  UT_ASSERT_EQUAL(C(1,1), 50.);
}



void
GEMMTest::testRectTransposedColMajor()
{
  Matrix<double> A(3,2), B(3,2), C(2,2);
  A(0,0) = 1;  A(0,1) = 2;
  A(1,0) = 3;  A(1,1) = 4;
  A(2,0) = 5;  A(2,1) = 6;
  B(0,0) = 7;  B(0,1) = 8;
  B(1,0) = 9;  B(1,1) = 10;
  B(2,0) = 11; B(2,1) = 12;
  C(0,0) = 0;  C(0,1) = 0;
  C(1,0) = 0;  C(1,1) = 0;

  Blas::gemm(1., A.t(), B, 0., C);

  UT_ASSERT_EQUAL(C(0,0), 89.);
  UT_ASSERT_EQUAL(C(0,1), 98.);
  UT_ASSERT_EQUAL(C(1,0), 116.);
  UT_ASSERT_EQUAL(C(1,1), 128.);
}


void
GEMMTest::testSquareColMajor()
{
  Matrix<double> A(2,2, false), B(2,2, false), C(2,2, false);
  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;
  B(0,0) = 5; B(0,1) = 6;
  B(1,0) = 7; B(1,1) = 8;
  C(0,0) = 0; C(0,1) = 0;
  C(1,0) = 0; C(1,1) = 0;

  Blas::gemm(1., A, B, 0., C);

  UT_ASSERT_EQUAL(C(0,0), 19.);
  UT_ASSERT_EQUAL(C(0,1), 22.);
  UT_ASSERT_EQUAL(C(1,0), 43.);
  UT_ASSERT_EQUAL(C(1,1), 50.);
}



UnitTest::TestSuite *
GEMMTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::gemm()");

  s->addTest(new UnitTest::TestCaller<GEMMTest>(
               "Blas::gemm(double[m,n], double[n,k], double[m,k]) (row-major)",
               &GEMMTest::testRectRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMMTest>(
               "Blas::gemm(double[m,n], double[n,k], double[m,k]) (col-major)",
               &GEMMTest::testRectColMajor));

  s->addTest(new UnitTest::TestCaller<GEMMTest>(
               "Blas::gemm(double[n,m]::t(), double[n,m], double[m,m]) (row-major)",
               &GEMMTest::testRectTransposedRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMMTest>(
               "Blas::gemm(double[m,m], double[m,m], double[m,m]) (row-major)",
               &GEMMTest::testSquareRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMMTest>(
               "Blas::gemm(double[n,m]::t(), double[n,m], double[m,m]) (col-major)",
               &GEMMTest::testRectTransposedColMajor));

  s->addTest(new UnitTest::TestCaller<GEMMTest>(
               "Blas::gemm(double[m,m], double[m,m], double[m,m]) (col-major)",
               &GEMMTest::testSquareColMajor));

  return s;
}
