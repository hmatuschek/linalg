#include "geqrftest.hh"
#include <lapack/geqrf.hh>
#include <lapack/ormqr.hh>

using namespace Linalg;


void
GEQRFTest::setUp()
{
  // Allocate some space;
  tau = Vector<double>(3);

  static double a_data[9] = { 12, -51,   4,
                               6, 167, -68,
                              -4,  24, -41 };

  A = Matrix<double>::fromData(a_data, 3, 3, 3, 1, 0).copy();
}



void
GEQRFTest::testRealNN()
{
  Vector<double> tmp(this->A.rows());
  Matrix<double> A = this->A.copy();
  Lapack::geqrf(A, tau);

  {
    Matrix<double> B = this->A.copy();
    Lapack::ormqr(A, tau, B, tmp, true, true);
    for (size_t i=0; i<3; i++) {
      for (size_t j=i; j<3; j++) {
        UT_ASSERT_NEAR(B(i,j), A(i,j));
      }
    }
  }

  {
    Matrix<double> B = this->A.copy().t();
    Lapack::ormqr(A, tau, B, tmp, false, false);
    for (size_t i=0; i<3; i++) {
      for (size_t j=i; j<3; j++) {
        UT_ASSERT_NEAR(B(j,i), A(i,j));
      }
    }
  }
}



UnitTest::TestSuite *
GEQRFTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for geqrf :");

  s->addTest(new UnitTest::TestCaller<GEQRFTest>(
               "Lapack::geqrf(double[n,n])", &GEQRFTest::testRealNN));

  return s;
}
