#include "dottest.hh"

#include "vector.hh"
#include "matrix.hh"
#include "blas/dot.hh"
#include "operators.hh"

using namespace Linalg;


void
DOTTest::testIncrFloat()
{
  Matrix<float> A = Matrix<float>::empty(5,2, true);

  A(0,0) = 1; A(1,0) = 2; A(2,0) = 3; A(3,0) = 4; A(4,0) = 5;
  A(0,1) = 6; A(1,1) = 7; A(2,1) = 8; A(3,1) = 9; A(4,1) = 10;

  UT_ASSERT_EQUAL(Blas::dot(A.col(0), A.col(1)), 130.f);
  UT_ASSERT_EQUAL(Blas::dot(A.col(0), A.col(1)), Blas::dot(A.col(1),A.col(0)));
}


void
DOTTest::testIncrDouble()
{
  Matrix<double> A = Matrix<double>::empty(5,2, true);

  A(0,0) = 1; A(1,0) = 2; A(2,0) = 3; A(3,0) = 4; A(4,0) = 5;
  A(0,1) = 6; A(1,1) = 7; A(2,1) = 8; A(3,1) = 9; A(4,1) = 10;

  UT_ASSERT_EQUAL(Blas::dot(A.col(0), A.col(1)), 130.);
  UT_ASSERT_EQUAL(Blas::dot(A.col(0), A.col(1)), Blas::dot(A.col(1),A.col(0)));
}


void
DOTTest::testDenseFloat()
{
  Matrix<float> A = Matrix<float>::empty(5,2, false);

  A(0,0) = 1; A(1,0) = 2; A(2,0) = 3; A(3,0) = 4; A(4,0) = 5;
  A(0,1) = 6; A(1,1) = 7; A(2,1) = 8; A(3,1) = 9; A(4,1) = 10;

  UT_ASSERT_EQUAL(Blas::dot(A.col(0), A.col(1)), 130.f);
  UT_ASSERT_EQUAL(Blas::dot(A.col(0), A.col(1)), Blas::dot(A.col(1),A.col(0)));
}


void
DOTTest::testDenseDouble()
{
  Matrix<double> A = Matrix<double>::empty(5,2, false);

  A(0,0) = 1; A(1,0) = 2; A(2,0) = 3; A(3,0) = 4; A(4,0) = 5;
  A(0,1) = 6; A(1,1) = 7; A(2,1) = 8; A(3,1) = 9; A(4,1) = 10;

  UT_ASSERT_EQUAL(Blas::dot(A.col(0), A.col(1)), 130.);
  UT_ASSERT_EQUAL(Blas::dot(A.col(0), A.col(1)), Blas::dot(A.col(1),A.col(0)));
}


void DOTTest::testHugeIncr()
{
  Matrix<float> A = Matrix<float>::empty(100*1024,2, true);
  Vector<float> x = A.col(0), y = A.col(1);

  for(size_t i=0; i<100; i++)
    Blas::dot(x, y);
}


void DOTTest::testHugeDense()
{
  Matrix<float> A = Matrix<float>::empty(100*1024,2, false);
  Vector<float> x = A.col(0), y = A.col(1);

  for(size_t i=0; i<100; i++)
    Blas::dot(x,y);
}


UnitTest::TestSuite *
DOTTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::dot()");

  s->addTest(new UnitTest::TestCaller<DOTTest>(
               "Blas::dot(float[m], float[m]) (incr)", &DOTTest::testIncrFloat));

  s->addTest(new UnitTest::TestCaller<DOTTest>(
               "Blas::dot(double[m], double[m]) (incr)", &DOTTest::testIncrDouble));

  s->addTest(new UnitTest::TestCaller<DOTTest>(
               "Blas::dot(float[m], float[m]) (dense)", &DOTTest::testDenseFloat));

  s->addTest(new UnitTest::TestCaller<DOTTest>(
               "Blas::dot(double[m], double[m]) (dense)", &DOTTest::testDenseDouble));

  s->addTest(new UnitTest::TestCaller<DOTTest>(
               "Blas::dot(double[m], double[m]) (huge-incr)", &DOTTest::testHugeIncr));

  s->addTest(new UnitTest::TestCaller<DOTTest>(
               "Blas::dot(double[m], double[m]) (huge-dense)", &DOTTest::testHugeDense));

  return s;
}
