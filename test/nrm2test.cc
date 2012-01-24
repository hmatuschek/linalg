#include "nrm2test.hh"

#include "vector.hh"
#include "matrix.hh"
#include "blas/nrm2.hh"

#include <cmath>

using namespace Linalg;


void
NRM2Test::testVector()
{
  // Simple direct vector test:
  Vector<double> vec(5);
  vec(0) = 1; vec(1) = 2; vec(2) = 3; vec(3) = 4; vec(4) = 5;

  // Test full vector:
  UT_ASSERT_EQUAL(Blas::nrm2(vec), std::sqrt(55));

  // test tail(1):
  UT_ASSERT_EQUAL(Blas::nrm2(vec.sub(1,4)), std::sqrt(54));

  // test head(4):
  UT_ASSERT_EQUAL(Blas::nrm2(vec.sub(0,4)), std::sqrt(30));
}


void
NRM2Test::testMatrixColumnRowMajor()
{
  Matrix<double> A(5,5);

  A(0,0) =  1, A(0, 1) =  2; A(0, 2) =  3; A(0, 3) =  4; A(0, 4) =  5;
  A(1,0) =  6, A(1, 1) =  7; A(1, 2) =  8; A(1, 3) =  9; A(1, 4) = 10;
  A(2,0) = 11, A(2, 1) = 12; A(2, 2) = 13; A(2, 3) = 14; A(2, 4) = 15;
  A(3,0) = 16, A(3, 1) = 17; A(3, 2) = 18; A(3, 3) = 19; A(3, 4) = 20;
  A(4,0) = 21, A(4, 1) = 22; A(4, 2) = 23; A(4, 3) = 24; A(4, 4) = 25;

  UT_ASSERT_EQUAL(Blas::nrm2(A.col(0)), std::sqrt(855));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(1)), std::sqrt(970));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(2)), std::sqrt(1095));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(3)), std::sqrt(1230));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(4)), std::sqrt(1375));

  UT_ASSERT_EQUAL(Blas::nrm2(A.col(0).sub(1,4)), std::sqrt(854));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(1).sub(1,4)), std::sqrt(966));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(2).sub(1,4)), std::sqrt(1086));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(3).sub(1,4)), std::sqrt(1214));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(4).sub(1,4)), std::sqrt(1350));

  UT_ASSERT_EQUAL(Blas::nrm2(A.col(0).sub(0,4)), std::sqrt(414));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(1).sub(0,4)), std::sqrt(486));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(2).sub(0,4)), std::sqrt(566));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(3).sub(0,4)), std::sqrt(654));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(4).sub(0,4)), std::sqrt(750));
}


void
NRM2Test::testMatrixRowRowMajor()
{
  Matrix<double> A(5,5);

  A(0,0) =  1, A(0, 1) =  2; A(0, 2) =  3; A(0, 3) =  4; A(0, 4) =  5;
  A(1,0) =  6, A(1, 1) =  7; A(1, 2) =  8; A(1, 3) =  9; A(1, 4) = 10;
  A(2,0) = 11, A(2, 1) = 12; A(2, 2) = 13; A(2, 3) = 14; A(2, 4) = 15;
  A(3,0) = 16, A(3, 1) = 17; A(3, 2) = 18; A(3, 3) = 19; A(3, 4) = 20;
  A(4,0) = 21, A(4, 1) = 22; A(4, 2) = 23; A(4, 3) = 24; A(4, 4) = 25;

  UT_ASSERT_EQUAL(Blas::nrm2(A.row(0)), std::sqrt(55));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(1)), std::sqrt(330));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(2)), std::sqrt(855));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(3)), std::sqrt(1630));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(4)), std::sqrt(2655));

  UT_ASSERT_EQUAL(Blas::nrm2(A.row(0).sub(1,4)), std::sqrt(54));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(1).sub(1,4)), std::sqrt(294));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(2).sub(1,4)), std::sqrt(734));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(3).sub(1,4)), std::sqrt(1374));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(4).sub(1,4)), std::sqrt(2214));

  UT_ASSERT_EQUAL(Blas::nrm2(A.row(0).sub(0,4)), std::sqrt(30));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(1).sub(0,4)), std::sqrt(230));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(2).sub(0,4)), std::sqrt(630));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(3).sub(0,4)), std::sqrt(1230));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(4).sub(0,4)), std::sqrt(2030));
}


void
NRM2Test::testMatrixColumnColMajor()
{
  Matrix<double> A(5,5, false);

  A(0,0) =  1, A(0, 1) =  2; A(0, 2) =  3; A(0, 3) =  4; A(0, 4) =  5;
  A(1,0) =  6, A(1, 1) =  7; A(1, 2) =  8; A(1, 3) =  9; A(1, 4) = 10;
  A(2,0) = 11, A(2, 1) = 12; A(2, 2) = 13; A(2, 3) = 14; A(2, 4) = 15;
  A(3,0) = 16, A(3, 1) = 17; A(3, 2) = 18; A(3, 3) = 19; A(3, 4) = 20;
  A(4,0) = 21, A(4, 1) = 22; A(4, 2) = 23; A(4, 3) = 24; A(4, 4) = 25;

  UT_ASSERT_EQUAL(Blas::nrm2(A.col(0)), std::sqrt(855));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(1)), std::sqrt(970));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(2)), std::sqrt(1095));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(3)), std::sqrt(1230));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(4)), std::sqrt(1375));

  UT_ASSERT_EQUAL(Blas::nrm2(A.col(0).sub(1,4)), std::sqrt(854));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(1).sub(1,4)), std::sqrt(966));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(2).sub(1,4)), std::sqrt(1086));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(3).sub(1,4)), std::sqrt(1214));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(4).sub(1,4)), std::sqrt(1350));

  UT_ASSERT_EQUAL(Blas::nrm2(A.col(0).sub(0,4)), std::sqrt(414));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(1).sub(0,4)), std::sqrt(486));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(2).sub(0,4)), std::sqrt(566));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(3).sub(0,4)), std::sqrt(654));
  UT_ASSERT_EQUAL(Blas::nrm2(A.col(4).sub(0,4)), std::sqrt(750));
}


void
NRM2Test::testMatrixRowColMajor()
{
  Matrix<double> A(5,5, false);

  A(0,0) =  1, A(0, 1) =  2; A(0, 2) =  3; A(0, 3) =  4; A(0, 4) =  5;
  A(1,0) =  6, A(1, 1) =  7; A(1, 2) =  8; A(1, 3) =  9; A(1, 4) = 10;
  A(2,0) = 11, A(2, 1) = 12; A(2, 2) = 13; A(2, 3) = 14; A(2, 4) = 15;
  A(3,0) = 16, A(3, 1) = 17; A(3, 2) = 18; A(3, 3) = 19; A(3, 4) = 20;
  A(4,0) = 21, A(4, 1) = 22; A(4, 2) = 23; A(4, 3) = 24; A(4, 4) = 25;

  UT_ASSERT_EQUAL(Blas::nrm2(A.row(0)), std::sqrt(55));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(1)), std::sqrt(330));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(2)), std::sqrt(855));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(3)), std::sqrt(1630));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(4)), std::sqrt(2655));

  UT_ASSERT_EQUAL(Blas::nrm2(A.row(0).sub(1,4)), std::sqrt(54));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(1).sub(1,4)), std::sqrt(294));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(2).sub(1,4)), std::sqrt(734));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(3).sub(1,4)), std::sqrt(1374));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(4).sub(1,4)), std::sqrt(2214));

  UT_ASSERT_EQUAL(Blas::nrm2(A.row(0).sub(0,4)), std::sqrt(30));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(1).sub(0,4)), std::sqrt(230));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(2).sub(0,4)), std::sqrt(630));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(3).sub(0,4)), std::sqrt(1230));
  UT_ASSERT_EQUAL(Blas::nrm2(A.row(4).sub(0,4)), std::sqrt(2030));
}



/*
 * Construct TestSuite:
 */
UnitTest::TestSuite *
NRM2Test::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::nrm2()");

  s->addTest(new UnitTest::TestCaller<NRM2Test>(
               "Blas::nrm2(double[m])", &NRM2Test::testVector));

  s->addTest(new UnitTest::TestCaller<NRM2Test>(
               "Blas::nrm2(double[m,m]::col(i) (row-major))", &NRM2Test::testMatrixColumnRowMajor));

  s->addTest(new UnitTest::TestCaller<NRM2Test>(
               "Blas::nrm2(double[m,m]::row(i) (row-major))", &NRM2Test::testMatrixRowRowMajor));

  s->addTest(new UnitTest::TestCaller<NRM2Test>(
               "Blas::nrm2(double[m,m]::col(i) (col-major))", &NRM2Test::testMatrixColumnColMajor));

  s->addTest(new UnitTest::TestCaller<NRM2Test>(
               "Blas::nrm2(double[m,m]::row(i) (col-major))", &NRM2Test::testMatrixRowColMajor));

  return s;
}
