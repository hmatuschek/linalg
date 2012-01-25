#include "matrixtest.hh"
#include "matrix.hh"

using namespace Linalg;


void
MatrixTest::testRectRowToColMajor()
{
  Matrix<double> A(2, 3);
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;

  Matrix<double> B(A.colMajor());

  for (size_t i=0; i<2; i++)
  {
    for (size_t j=0; j<3; j++)
    {
      UT_ASSERT_EQUAL(A(i,j), B(i,j));
      UT_ASSERT_EQUAL(A.t()(j,i), B(i,j));
      UT_ASSERT_EQUAL(A(i,j), B.t()(j,i));
      UT_ASSERT_EQUAL(A.t()(j,i), B.t()(j,i));
    }
  }
}


void
MatrixTest::testSquareRowToColMajor()
{
  Matrix<double> A(2, 2);
  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;

  Matrix<double> B(A.colMajor());

  for (size_t i=0; i<2; i++)
  {
    for (size_t j=0; j<2; j++)
    {
      UT_ASSERT_EQUAL(A(i,j), B(i,j));
      UT_ASSERT_EQUAL(A.t()(j,i), B(i,j));
      UT_ASSERT_EQUAL(A(i,j), B.t()(j,i));
      UT_ASSERT_EQUAL(A.t()(j,i), B.t()(j,i));
    }
  }
}



void
MatrixTest::testRectColToRowMajor()
{
  Matrix<double> A(2, 3, false);
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;

  Matrix<double> B(A.rowMajor());

  for (size_t i=0; i<2; i++)
  {
    for (size_t j=0; j<3; j++)
    {
      UT_ASSERT_EQUAL(A(i,j), B(i,j));
      UT_ASSERT_EQUAL(A.t()(j,i), B(i,j));
      UT_ASSERT_EQUAL(A(i,j), B.t()(j,i));
      UT_ASSERT_EQUAL(A.t()(j,i), B.t()(j,i));
    }
  }
}


void
MatrixTest::testSquareColToRowMajor()
{
  Matrix<double> A(2, 2, false);
  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;

  Matrix<double> B(A.rowMajor());

  for (size_t i=0; i<2; i++)
  {
    for (size_t j=0; j<2; j++)
    {
      UT_ASSERT_EQUAL(A(i,j), B(i,j));
      UT_ASSERT_EQUAL(A.t()(j,i), B(i,j));
      UT_ASSERT_EQUAL(A(i,j), B.t()(j,i));
      UT_ASSERT_EQUAL(A.t()(j,i), B.t()(j,i));
    }
  }
}



void
MatrixTest::testFromData()
{
  double A_data[6] = {1., 2., 3., 4., 5., 6.};
  Matrix<double> A = matrixFromData<>(A_data, 2, 3);

  double B_data[6] = {1., 4., 2., 5., 3., 6.};
  Matrix<double> B = matrixFromData<>(B_data, 2, 3, 0, 0, false, false, false);

  for (size_t i=0; i<A.rows(); i++) {
    for (size_t j=0; j<A.cols(); j++) {
      UT_ASSERT_EQUAL(A(i,j), B(i,j));
      UT_ASSERT_EQUAL(A.colMajor()(i,j), B(i,j));
      UT_ASSERT_EQUAL(A(i,j), B.rowMajor()(i,j));
      UT_ASSERT_EQUAL(A.colMajor()(i,j), B.rowMajor()(i,j));
    }
  }
}



UnitTest::TestSuite *
MatrixTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Matrix class");

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,m]::colMajor() (row-major)", &MatrixTest::testSquareRowToColMajor));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,n]::colMajor() (row-major)", &MatrixTest::testRectRowToColMajor));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,m]::rowMajor() (col-major)", &MatrixTest::testSquareColToRowMajor));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,n]::rowMajor() (col-major)", &MatrixTest::testRectColToRowMajor));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,n] internal data layout.", &MatrixTest::testFromData));

  return s;
}
