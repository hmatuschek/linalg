#include "trimatrixtest.hh"

#include "trimatrix.hh"

using namespace Linalg;



void
TriMatrixTest::testTranspose()
{
  Matrix<double> A(3,3);
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;
  A(2,0) = 7; A(2,1) = 8; A(2,2) = 9;

  UT_ASSERT(! triu(A).t().isUpper());
  UT_ASSERT(tril(A).t().isUpper());
  UT_ASSERT(triu(A.t()).isUpper());
  UT_ASSERT(! tril(A.t()).isUpper());

  for (size_t i=0; i<A.rows(); i++) {
    for (size_t j=i; j<A.cols(); j++) {
      UT_ASSERT_EQUAL(triu(A).t()(j,i), A(i,j));
      UT_ASSERT_EQUAL(triu(A.t())(i,j), A(j,i));
    }
  }
}



void
TriMatrixTest::testRowToColMajor()
{
  Matrix<double> A(3,3);
  Matrix<double> B(A.colMajor());
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;
  A(2,0) = 7; A(2,1) = 8; A(2,2) = 9;

  TriMatrix<double> Au(triu(A));
  TriMatrix<double> Bu(triu(B));
  for(size_t i=0; i<A.rows(); i++)
  {
    for (size_t j=i; j<A.cols(); j++)
    {
      UT_ASSERT_EQUAL(Au(i,j), A(i,j));
      UT_ASSERT_EQUAL(Bu(i,j), A(i,j));
    }
  }

  TriMatrix<double> Al(tril(A));
  TriMatrix<double> Bl(tril(B));
  for(size_t j=0; j<A.cols(); j++)
  {
    for (size_t i=j; i<A.rows(); i++)
    {
      UT_ASSERT_EQUAL(Al(i,j), A(i,j));
      UT_ASSERT_EQUAL(Bl(i,j), A(i,j));
    }
  }
}



void
TriMatrixTest::testColToRowMajor()
{
  Matrix<double> A(3,3, false);
  Matrix<double> B(A.rowMajor());
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;
  A(2,0) = 7; A(2,1) = 8; A(2,2) = 9;

  TriMatrix<double> Au(triu(A));
  TriMatrix<double> Bu(triu(B));
  for(size_t i=0; i<A.rows(); i++)
  {
    for (size_t j=i; j<A.cols(); j++)
    {
      UT_ASSERT_EQUAL(Au(i,j), A(i,j));
      UT_ASSERT_EQUAL(Bu(i,j), A(i,j));
    }
  }

  TriMatrix<double> Al(tril(A));
  TriMatrix<double> Bl(tril(B));
  for(size_t j=0; j<A.cols(); j++)
  {
    for (size_t i=j; i<A.rows(); i++)
    {
      UT_ASSERT_EQUAL(Al(i,j), A(i,j));
      UT_ASSERT_EQUAL(Bl(i,j), A(i,j));
    }
  }
}



UnitTest::TestSuite *
TriMatrixTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for TriMatrix class");

  s->addTest(new UnitTest::TestCaller<TriMatrixTest>(
               "triu(double[m,m])::t()", &TriMatrixTest::testTranspose));

  s->addTest(new UnitTest::TestCaller<TriMatrixTest>(
               "triu(double[m,m]::colMajor()) (row-major)", &TriMatrixTest::testRowToColMajor));

  s->addTest(new UnitTest::TestCaller<TriMatrixTest>(
               "triu(double[m,m]::rowMajor()) (col-major)", &TriMatrixTest::testColToRowMajor));

  return s;
}
