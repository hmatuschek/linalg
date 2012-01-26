#include "trmvtest.hh"

#include "trimatrix.hh"
#include "vector.hh"
#include "blas/trmv.hh"

using namespace Linalg;


void
TRMVTest::setUp()
{
  this->Arow = Matrix<double>::empty(3,3);
  this->Acol = Matrix<double>::empty(3,3, false);
  this->x    = Vector<double>::empty(3);

  Matrix<double> A(3,3);
  A(0,0)=1; A(0,1)=2; A(0,2)=3;
  A(1,0)=4; A(1,1)=5; A(1,2)=6;
  A(2,0)=7; A(2,1)=8; A(2,2)=9;

  Vector<double> x(3);
  x(0) = 10; x(1) = 11; x(2) = 12;

  this->Arow.vals() = A;
  this->Acol.vals() = A;
  this->x.vals()    = x;
}


void
TRMVTest::testUpperRowMajor()
{
  Vector<double> x(3); x.vals() = this->x;

  Blas::trmv(triu(this->Arow), x);
  UT_ASSERT_EQUAL(x(0), 1*10+2*11+3*12.);
  UT_ASSERT_EQUAL(x(1),      5*11+6*12.);
  UT_ASSERT_EQUAL(x(2),           9*12.);
}


void
TRMVTest::testUpperTransRowMajor()
{
  Vector<double> x(3); x.vals() = this->x;

  Blas::trmv(triu(this->Arow).t(), x);
  UT_ASSERT_EQUAL(x(0), 1*10+0*11+0*12.);
  UT_ASSERT_EQUAL(x(1), 2*10+5*11+0*12.);
  UT_ASSERT_EQUAL(x(2), 3*10+6*11+9*12.);
}


void
TRMVTest::testTransUpperRowMajor()
{
  Vector<double> x(3); x.vals() = this->x;

  Blas::trmv(triu(this->Arow.t()), x);
  UT_ASSERT_EQUAL(x(0), 1*10+4*11+7*12.);
  UT_ASSERT_EQUAL(x(1),      5*11+8*12.);
  UT_ASSERT_EQUAL(x(2),           9*12.);
}


void
TRMVTest::testUpperColMajor()
{
  Vector<double> x(3); x.vals() = this->x;

  Blas::trmv(triu(this->Acol), x);
  UT_ASSERT_EQUAL(x(0), 1*10+2*11+3*12.);
  UT_ASSERT_EQUAL(x(1),      5*11+6*12.);
  UT_ASSERT_EQUAL(x(2),           9*12.);
}


void
TRMVTest::testUpperTransColMajor()
{
  Vector<double> x(3); x.vals() = this->x;

  Blas::trmv(triu(this->Acol).t(), x);
  UT_ASSERT_EQUAL(x(0), 1*10+0*11+0*12.);
  UT_ASSERT_EQUAL(x(1), 2*10+5*11+0*12.);
  UT_ASSERT_EQUAL(x(2), 3*10+6*11+9*12.);
}


void
TRMVTest::testTransUpperColMajor()
{
  Vector<double> x(3); x.vals() = this->x;

  Blas::trmv(triu(this->Acol.t()), x);
  UT_ASSERT_EQUAL(x(0), 1*10+4*11+7*12.);
  UT_ASSERT_EQUAL(x(1),      5*11+8*12.);
  UT_ASSERT_EQUAL(x(2),           9*12.);
}



UnitTest::TestSuite *
TRMVTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Blas::trmv()");

  s->addTest(new UnitTest::TestCaller<TRMVTest>(
               "Blas::trmv(triu(double[m,m]), double[m]) (row-major)",
               &TRMVTest::testUpperRowMajor));

  s->addTest(new UnitTest::TestCaller<TRMVTest>(
               "Blas::trmv(triu(double[m,m])::t(), double[m]) (row-major)",
               &TRMVTest::testUpperTransRowMajor));

  s->addTest(new UnitTest::TestCaller<TRMVTest>(
               "Blas::trmv(triu(double[m,m]::t()), double[m]) (row-major)",
               &TRMVTest::testTransUpperRowMajor));

  s->addTest(new UnitTest::TestCaller<TRMVTest>(
               "Blas::trmv(triu(double[m,m]), double[m]) (col-major)",
               &TRMVTest::testUpperColMajor));

  s->addTest(new UnitTest::TestCaller<TRMVTest>(
               "Blas::trmv(triu(double[m,m])::t(), double[m]) (col-major)",
               &TRMVTest::testUpperColMajor));

  s->addTest(new UnitTest::TestCaller<TRMVTest>(
               "Blas::trmv(triu(double[m,m]::t()), double[m]) (col-major)",
               &TRMVTest::testUpperColMajor));

  return s;
}
