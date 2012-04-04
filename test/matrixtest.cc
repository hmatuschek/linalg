#include "matrixtest.hh"
#include "matrix.hh"

using namespace Linalg;


void
MatrixTest::testFromData()
{
  double A_data[6] = {1., 2., 3., 4., 5., 6.};
  Matrix<double> A = Matrix<double>::fromData(A_data, 2, 3, 3, 1, 0);

  double B_data[6] = {1., 4., 2., 5., 3., 6.};
  Matrix<double> B = Matrix<double>::fromData(B_data, 2, 3, 1, 2, 0);

  UT_ASSERT_EQUAL((int)A.rows(), 2); UT_ASSERT_EQUAL((int)A.cols(), 3);
  UT_ASSERT_EQUAL((int)B.rows(), 2); UT_ASSERT_EQUAL((int)B.cols(), 3);
  UT_ASSERT_EQUAL((int)A.strides(0), 3); UT_ASSERT_EQUAL((int)A.strides(1), 1);
  UT_ASSERT_EQUAL((int)B.strides(0), 1); UT_ASSERT_EQUAL((int)B.strides(1), 2);

  for (size_t i=0; i<A.rows(); i++) {
    for (size_t j=0; j<A.cols(); j++) {
      UT_ASSERT_EQUAL(A(i,j), B(i,j));
      UT_ASSERT_EQUAL(A.t()(j,i), B(i,j));
    }
  }
}



void
MatrixTest::testValueRef()
{
  Matrix<double> A = Matrix<double>::rand(3,2);
  Matrix<double> B(3,2);

  // Assign values:
  B.values() = A;

  // compare values:
  for (size_t i=0; i<A.rows(); i++) {
    for (size_t j=0; j<A.cols(); j++) {
      UT_ASSERT_EQUAL(B(i,j), A(i,j));
    }
  }
}


void
MatrixTest::testRowSelRowMajor()
{
  double data[9] = {1,2,3,4,5,6,7,8,9};
  Matrix<double> A = Matrix<double>::fromData(data, 3, 3, 3, 1, 0);
  Vector<double> a;

  a = A.row(0);
  UT_ASSERT_EQUAL(a(0), 1.); UT_ASSERT_EQUAL(a(1), 2.); UT_ASSERT_EQUAL(a(2), 3.);
  a = A.row(1);
  UT_ASSERT_EQUAL(a(0), 4.); UT_ASSERT_EQUAL(a(1), 5.); UT_ASSERT_EQUAL(a(2), 6.);
  a = A.row(2);
  UT_ASSERT_EQUAL(a(0), 7.); UT_ASSERT_EQUAL(a(1), 8.); UT_ASSERT_EQUAL(a(2), 9.);
}


void
MatrixTest::testRowSelColMajor()
{
  double data[9] = {1,2,3,4,5,6,7,8,9};
  Matrix<double> A = Matrix<double>::fromData(data, 3, 3, 1, 3, 0);
  Vector<double> a;

  a = A.row(0);
  UT_ASSERT_EQUAL(a(0), 1.); UT_ASSERT_EQUAL(a(1), 4.); UT_ASSERT_EQUAL(a(2), 7.);
  a = A.row(1);
  UT_ASSERT_EQUAL(a(0), 2.); UT_ASSERT_EQUAL(a(1), 5.); UT_ASSERT_EQUAL(a(2), 8.);
  a = A.row(2);
  UT_ASSERT_EQUAL(a(0), 3.); UT_ASSERT_EQUAL(a(1), 6.); UT_ASSERT_EQUAL(a(2), 9.);
}


void
MatrixTest::testColSelRowMajor()
{
  double data[9] = {1,2,3,4,5,6,7,8,9};
  Matrix<double> A = Matrix<double>::fromData(data, 3, 3, 3, 1, 0);
  Vector<double> a;

  a = A.col(0);
  UT_ASSERT_EQUAL(a(0), 1.); UT_ASSERT_EQUAL(a(1), 4.); UT_ASSERT_EQUAL(a(2), 7.);
  a = A.col(1);
  UT_ASSERT_EQUAL(a(0), 2.); UT_ASSERT_EQUAL(a(1), 5.); UT_ASSERT_EQUAL(a(2), 8.);
  a = A.col(2);
  UT_ASSERT_EQUAL(a(0), 3.); UT_ASSERT_EQUAL(a(1), 6.); UT_ASSERT_EQUAL(a(2), 9.);
}


void
MatrixTest::testColSelColMajor()
{
  double data[9] = {1,2,3,4,5,6,7,8,9};
  Matrix<double> A = Matrix<double>::fromData(data, 3, 3, 1, 3, 0);
  Vector<double> a;

  a = A.col(0);
  UT_ASSERT_EQUAL(a(0), 1.); UT_ASSERT_EQUAL(a(1), 2.); UT_ASSERT_EQUAL(a(2), 3.);
  a = A.col(1);
  UT_ASSERT_EQUAL(a(0), 4.); UT_ASSERT_EQUAL(a(1), 5.); UT_ASSERT_EQUAL(a(2), 6.);
  a = A.col(2);
  UT_ASSERT_EQUAL(a(0), 7.); UT_ASSERT_EQUAL(a(1), 8.); UT_ASSERT_EQUAL(a(2), 9.);
}


void
MatrixTest::testSwap()
{
  double a_data[6] = {1,2, 3,4, 5,6};
  Matrix<double> A = Matrix<double>::fromData(a_data, 3,2, 2, 1, 0);
  double b_data[6] = {6,5,4, 3,2,1};
  Matrix<double> B = Matrix<double>::fromData(b_data, 2,3, 3, 1, 0);

  A.swap(B);
  UT_ASSERT_EQUAL((int)A.rows(), 2); UT_ASSERT_EQUAL((int)A.cols(), 3);
  UT_ASSERT_EQUAL((int)B.rows(), 3); UT_ASSERT_EQUAL((int)B.cols(), 2);

  UT_ASSERT_EQUAL(A(0,0), 6.); UT_ASSERT_EQUAL(A(0,1), 5.); UT_ASSERT_EQUAL(A(0,2), 4.);
  UT_ASSERT_EQUAL(A(1,0), 3.); UT_ASSERT_EQUAL(A(1,1), 2.); UT_ASSERT_EQUAL(A(1,2), 1.);

  UT_ASSERT_EQUAL(B(0,0), 1.); UT_ASSERT_EQUAL(B(0,1), 2.);
  UT_ASSERT_EQUAL(B(1,0), 3.); UT_ASSERT_EQUAL(B(1,1), 4.);
  UT_ASSERT_EQUAL(B(2,0), 5.); UT_ASSERT_EQUAL(B(2,1), 6.);
}


UnitTest::TestSuite *
MatrixTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Matrix class");

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,n] internal data layout", &MatrixTest::testFromData));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,n]::vals() value assignment", &MatrixTest::testValueRef));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,m]::row() (row-major)", &MatrixTest::testRowSelRowMajor));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,m]::row() (col-major)", &MatrixTest::testRowSelColMajor));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,m]::col() (row-major)", &MatrixTest::testColSelRowMajor));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,m]::col() (col-major)", &MatrixTest::testColSelColMajor));

  s->addTest(new UnitTest::TestCaller<MatrixTest>(
               "double[m,n]::swap(double[n,m])", &MatrixTest::testSwap));

  return s;
}
