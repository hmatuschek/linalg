#include "gemvtest.hh"
#include "matrix.hh"
#include "vector.hh"
#include "blas/gemv.hh"

using namespace Linalg;


void
GEMVTest::testSquareRowMajor()
{
  Matrix<double> A(3,3);
  Vector<double> x(3), y(3);

  A(0,0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1,0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  A(2,0) = 7; A(2, 1) = 8; A(2, 2) = 9;
  x(0) = 10; x(1) = 11; x(2) = 12;
  y(0) = 0; y(1) = 0; y(2) = 0;

  Blas::gemv(1, A, x, 0.0, y);

  UT_ASSERT_EQUAL(y(0), 68.);
  UT_ASSERT_EQUAL(y(1), 40.+55.+72.);
  UT_ASSERT_EQUAL(y(2), 70.+88.+108.);
}


void
GEMVTest::testSquareTransposedRowMajor()
{
  Matrix<double> A(3,3);
  Vector<double> x(3), y(3);

  A(0,0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1,0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  A(2,0) = 7; A(2, 1) = 8; A(2, 2) = 9;
  x(0) = 10; x(1) = 11; x(2) = 12;
  y(0) = 0; y(1) = 0; y(2) = 0;

  Blas::gemv(1, A.t(), x, 0.0, y);

  UT_ASSERT_EQUAL(y(0), 10. + 44. + 84.);
  UT_ASSERT_EQUAL(y(1), 20. + 55. + 96.);
  UT_ASSERT_EQUAL(y(2), 30. + 66. + 108.);
}


void
GEMVTest::testSquareColMajor()
{
  Matrix<double> A(3,3);
  Vector<double> x(3), y(3);

  A(0,0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1,0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  A(2,0) = 7; A(2, 1) = 8; A(2, 2) = 9;
  x(0) = 10; x(1) = 11; x(2) = 12;
  y(0) = 0; y(1) = 0; y(2) = 0;

  Blas::gemv(1, A, x, 0.0, y);

  UT_ASSERT_EQUAL(y(0), 68.);
  UT_ASSERT_EQUAL(y(1), 40.+55.+72.);
  UT_ASSERT_EQUAL(y(2), 70.+88.+108.);
}


void
GEMVTest::testSquareTransposedColMajor()
{
  Matrix<double> A(3,3);
  Vector<double> x(3), y(3);

  A(0,0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1,0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  A(2,0) = 7; A(2, 1) = 8; A(2, 2) = 9;
  x(0) = 10; x(1) = 11; x(2) = 12;
  y(0) = 0; y(1) = 0; y(2) = 0;

  Blas::gemv(1, A.t(), x, 0.0, y);

  UT_ASSERT_EQUAL(y(0), 10. + 44. + 84.);
  UT_ASSERT_EQUAL(y(1), 20. + 55. + 96.);
  UT_ASSERT_EQUAL(y(2), 30. + 66. + 108.);
}



void
GEMVTest::testRectRowMajor()
{
  Matrix<double> A(2,3);
  Vector<double> x(3), y(2);

  A(0,0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1,0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  x(0) = 10; x(1) = 11; x(2) = 12;
  y(0) = 0; y(1) = 0;

  Blas::gemv(1, A, x, 0.0, y);

  UT_ASSERT_EQUAL(y(0), 68.);
  UT_ASSERT_EQUAL(y(1), 167.);
}


void
GEMVTest::testRectTransposedRowMajor()
{
  Matrix<double> A(2,3);
  Vector<double> x(2), y(3);

  A(0,0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1,0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  x(0) = 10; x(1) = 11;
  y(0) = 0; y(1) = 0; y(2) = 0;

  Blas::gemv(1, A.t(), x, 0.0, y);

  UT_ASSERT_EQUAL(y(0), 54.);
  UT_ASSERT_EQUAL(y(1), 75.);
  UT_ASSERT_EQUAL(y(2), 96.);
}


void
GEMVTest::testRectColMajor()
{
  Matrix<double> A(2,3);
  Vector<double> x(3), y(2);

  A(0,0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1,0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  x(0) = 10; x(1) = 11; x(2) = 12;
  y(0) = 0; y(1) = 0;

  Blas::gemv(1, A, x, 0.0, y);

  UT_ASSERT_EQUAL(y(0), 68.);
  UT_ASSERT_EQUAL(y(1), 167.);
}


void
GEMVTest::testRectTransposedColMajor()
{
  Matrix<double> A(2,3);
  Vector<double> x(2), y(3);

  A(0,0) = 1; A(0, 1) = 2; A(0, 2) = 3;
  A(1,0) = 4; A(1, 1) = 5; A(1, 2) = 6;
  x(0) = 10; x(1) = 11;
  y(0) = 0; y(1) = 0; y(2) = 0;

  Blas::gemv(1, A.t(), x, 0.0, y);

  UT_ASSERT_EQUAL(y(0), 54.);
  UT_ASSERT_EQUAL(y(1), 75.);
  UT_ASSERT_EQUAL(y(2), 96.);
}



UnitTest::TestSuite *
GEMVTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::gemv()");

  s->addTest(new UnitTest::TestCaller<GEMVTest>(
               "Blas::gemv(double[m,m], double[m], double[m]) (row-major)",
               &GEMVTest::testSquareRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>(
               "Blas::gemv(double[m,n], double[n], double[m]) (row-major)",
               &GEMVTest::testRectRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>(
               "Blas::gemv(double[m,m]::t(), double[m], double[m]) (row-major)",
               &GEMVTest::testSquareTransposedRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>(
               "Blas::gemv(double[m,n]::t(), double[m], double[n]) (row-major)",
               &GEMVTest::testRectTransposedRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>(
               "Blas::gemv(double[m,m], double[m], double[m]) (col-major)",
               &GEMVTest::testSquareColMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>(
               "Blas::gemv(double[m,n], double[n], double[m]) (col-major)",
               &GEMVTest::testRectColMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>(
               "Blas::gemv(double[m,m]::t(), double[m], double[m]) (col-major)",
               &GEMVTest::testSquareTransposedColMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>(
               "Blas::gemv(double[m,n]::t(), double[m], double[n]) (col-major)",
               &GEMVTest::testRectTransposedColMajor));

  return s;

}
