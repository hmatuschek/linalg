#include "gemvtest.hh"
#include "matrix.hh"
#include "vector.hh"
#include "blas/gemv.hh"

using namespace Linalg;


void
GEMVTest::testRowMajor()
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
GEMVTest::testTransposedRowMajor()
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
GEMVTest::testColMajor()
{
  Matrix<double> A(3,3, false);
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
GEMVTest::testTransposedColMajor()
{
  Matrix<double> A(3,3, false);
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



UnitTest::TestSuite *
GEMVTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::gemv()");

  s->addTest(new UnitTest::TestCaller<GEMVTest>("Blas::gemv(Matrix<double>, Vector<double>, Vector<double>) (row-major)",
                                                &GEMVTest::testRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>("Blas::gemv(Matrix<double>::t(), Vector<double>, Vector<double>) (row-major)",
                                                &GEMVTest::testTransposedRowMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>("Blas::gemv(Matrix<double>, Vector<double>, Vector<double>) (col-major)",
                                                &GEMVTest::testColMajor));

  s->addTest(new UnitTest::TestCaller<GEMVTest>("Blas::gemv(Matrix<double>::t(), Vector<double>, Vector<double>) (col-major)",
                                                &GEMVTest::testTransposedColMajor));

  return s;

}
