#include "dottest.hh"

#include "vector.hh"
#include "blas/dot.hh"

using namespace Linalg;


void
DOTTest::testVectorVector()
{
  Vector<double> a(5), b(5);

  a(0) = 1; a(1) = 2; a(2) = 3; a(3) = 4; a(4) = 5;
  b(0) = 6; b(1) = 7; b(2) = 8; b(3) = 9; b(4) = 10;

  UT_ASSERT_EQUAL(Blas::dot(a,b), 130.);
  UT_ASSERT_EQUAL(Blas::dot(b,a), Blas::dot(a,b));
}



UnitTest::TestSuite *
DOTTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::dot()");

  s->addTest(new UnitTest::TestCaller<DOTTest>(
               "Blas::dot(double[m], double[m])", &DOTTest::testVectorVector));

  return s;
}
